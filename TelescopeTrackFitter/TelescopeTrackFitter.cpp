/*
 * TelescopeTrackFitter.cpp
 *
 *  Created on: Jun 22, 2017
 *      Author: cligtenb
 */

#include "TelescopeTrackFitter.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iterator>
#include <list>
#include <utility>

#include "TFile.h"
#include "TTree.h"
#include "TVirtualPad.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TH1.h"
#include "TRandom.h"

#include "../Alignment/Alignment.h"
#include "DetectorConfiguration.h"
#include "/user/cligtenb/rootmacros/getObjectFromFile.h"
//#include "../../rootmacros/getObjectFromFile.h"
#include "ResidualHistogrammer.h"
#include "transformHits.h"
#include "makeNoisyPixelMask.cpp"

using namespace std;

namespace {

	struct TrackState {
		TVector3 position;
		double angle[2];
		constexpr static double pixelSize=18.4E-3;//mm

		TrackState(TVector3 position, double xangle, double yangle) : position(position) {
			angle[0]=xangle;
			angle[1]=yangle;
		}

		void scatter(double scatterAngle) {
			for(int i=0; i<2; i++) {
				angle[i]+=gRandom->Gaus(0,scatterAngle);
			}
		}
		void scatterInVolume(double scatterAngle, double distance) {
			//scatter using PDG approach
			propagate(distance);
			for(int i=0; i<2; i++) {
				double z1 = gRandom->Gaus(0,scatterAngle), z2=gRandom->Gaus(0,scatterAngle);
				angle[i]+=z2;
				position[i]+=z1*distance/sqrt(12)+z2*distance/2;
			}
		}

		void propagate(double distance) {
			for(int i=0; i<2; i++) {
				position[i]+=distance*angle[i];
			}
			position[2]+=distance;
		}
		void propagateTo(double zcoordinate) {
			propagate(zcoordinate-position.z());
		}

		PositionHit getHit(unsigned char plane) {
			 return {position.x()+gRandom->Uniform(-pixelSize/2, pixelSize/2), position.y()+gRandom->Uniform(-pixelSize/2, pixelSize/2), position.z(), plane};
		}

	};

	std::vector<std::vector<PositionHit> > generateSpaceHits(const TelescopeConfiguration& detector, TVector3& simulatedQuadPosition) {
		TrackState track{ {5 + gRandom->Uniform(-1,1), 5 + gRandom->Uniform(-1,1),0}, gRandom->Uniform(-1E-4,1E-4),gRandom->Uniform(-1E-4,1E-4)};
		std::vector<std::vector<PositionHit> > hits(6);

		const double scatterPerPlane=0.108E-3;
		const double quadEntry=140, quadExit=240;

		//for reference: {0, 19.7, 38.9, 396.6, 415.8, 435.4 }, //plane positions as measured in oktober18 quad

		track.propagateTo(detector.planePosition[0]);
//		track.scatter(scatterPerPlane);
		hits[0].push_back( track.getHit(0) );
		track.propagateTo(detector.planePosition[1]);
		track.scatter(scatterPerPlane);
		hits[1].push_back( track.getHit(1) );
		track.propagateTo(detector.planePosition[2]);
		track.scatter(scatterPerPlane);
		hits[2].push_back( track.getHit(2) );


//		track.propagateTo(220);
//		track.scatter(0.24E-3);
//		track.propagateTo(quadEntry);
//		track.propagateTo(quadExit);

		track.scatterInVolume(0.069E-3, quadEntry-detector.planePosition[2]);

		track.scatter(0.048E-3);
		track.scatterInVolume(0.070E-3, (quadExit-quadEntry)/2);
		simulatedQuadPosition=track.getHit(0).position;
		track.scatterInVolume(0.070E-3, (quadExit-quadEntry)/2);
		track.scatter(0.048E-3);
		track.scatter(0.16E-3);//extra scatter

		track.scatterInVolume(0.088E-3, detector.planePosition[3]-quadExit);

		track.propagateTo(detector.planePosition[3]);
		track.scatter(scatterPerPlane);
		hits[3].push_back( track.getHit(3) );
		track.propagateTo(detector.planePosition[4]);
		track.scatter(scatterPerPlane);
		hits[4].push_back( track.getHit(4) );
		track.propagateTo(detector.planePosition[5]);
//		track.scatter(scatterPerPlane);
		hits[5].push_back( track.getHit(5) );


		//set errors
		for(auto& hv: hits) for(auto& h : hv) h.error=TVector3(detector.pixelsize/sqrt(12),detector.pixelsize/sqrt(12),1);

		return hits;
	}

	TVector3 averagePosition(const std::list<PositionHit>& hl){
		TVector3 sum;
		for(auto& h : hl) sum+=h.position;
		return sum*(1.0/hl.size());
	}

	HoughTransformer::HitCluster clusterHitsPerPlane(const HoughTransformer::HitCluster& hc, const TelescopeConfiguration& detector) {
		std::vector<std::list<PositionHit>> hitsPerPlane(detector.nPlanes);
		for(const auto& h : hc) {
			hitsPerPlane.at(h.chip).push_back(h);
		}
		HoughTransformer::HitCluster clustered;
		for(unsigned char i=0; i<hitsPerPlane.size(); i++) {
			if(hitsPerPlane[i].empty()) continue;
			if(hitsPerPlane[i].size()==1) {
				clustered.add(hitsPerPlane[i].front());
			} else {
				auto av=averagePosition(hitsPerPlane[i]);
				clustered.add(PositionHit{av.x(), av.y(), av.z(), i});
			}
		}

		//set errors for the moment to the same
		for(auto& h : clustered) h.error=TVector3(detector.pixelsize/sqrt(12),detector.pixelsize/sqrt(12),1);

		return clustered;
	}

}//end namespace


TelescopeTrackFitter::TelescopeTrackFitter(std::string inputfile, const TelescopeConfiguration& detector) :
	detector(detector),
	houghTransform(detector.xmin(), detector.xmax(), detector.ymin(), detector.ymax(), 100/*xbins*/, 200 /*ybins*/ ), //DEBUG was 200 500
	binnedClustering(detector.xmin(), detector.xmax(), detector.ymin(), detector.ymax(),10,5),
	residualHistograms(nullptr),
	hitsCentre(detector.nPlanes, detector.getCentre() ), //initialise to regular centre of sensor
	averageResidualFromSum(detector.nPlanes, {0,0} ), //initialise to zero
	rotationZFromSum(detector.nPlanes),
	slope1FromSum(0), slope2FromSum(0)
{

	//get tree from file
	file=openFile(inputfile);
	hitTable=getObjectFromFile<TTree>("Hits", file);


	//setup tree for reading
	hitTable->SetBranchAddress("mimosa", &mimosaHit);
	nEvents=hitTable->GetEntries();
	//    unsigned short triggerNumberBegin, triggerNumberEnd;
	hitTable->SetBranchAddress("triggerNumberBegin", &triggerNumberBegin);
	hitTable->SetBranchAddress("triggerNumberEnd", &triggerNumberEnd);
	hitTable->SetBranchAddress("timestamp", &timestamp);

}

TelescopeTrackFitter::~TelescopeTrackFitter() {
	file->Close();
	residualHistograms=nullptr; //make sure to close residualHistograms before track, because they use the same file
	trackHistograms=nullptr; //todo: combine both in a proper way;
	simulatedQuadxResidual.release();//cannot delete once saved to a file?
	simulatedQuadyResidual.release();//cannot delete once saved to a file?
}

int TelescopeTrackFitter::makeMask(double ntimesThreshold) {
	//first mask noisy pixels
	for(int iplane=0; iplane<detector.nPlanes; iplane++) {
		mask.emplace_back( makeNoisyPixelMask(hitTable,iplane,ntimesThreshold,{detector.pixelColumns, detector.pixelRows} ) );
	}
	return 0; //should return number of pixels masked
}

bool TelescopeTrackFitter::passEvent(const std::vector<std::vector<PositionHit> >& spaceHit) const {
//	if(triggerNumberBegin==triggerNumberEnd) return false;
	bool planeIsHit[6]={0};
	int nPlanesHit=0;
	int nTotalHits=0;
	for(int i=0; i<detector.nPlanes; i++) {
		int size=spaceHit[i].size();
		if(size) {
			nPlanesHit++;
			planeIsHit[i]=true;
		}
		nTotalHits+=size;
	}
	const int nMinPlanesHitPerPart=2;
	return (planeIsHit[0]+planeIsHit[1]+planeIsHit[2])>=nMinPlanesHitPerPart and (planeIsHit[3]+planeIsHit[4]+planeIsHit[5])>=nMinPlanesHitPerPart;
}

std::vector<std::vector<PositionHit> > TelescopeTrackFitter::getSpaceHits() {
	std::vector<std::vector<PositionHit> > spaceHit;
	//apply mask
	if (mask.empty()) {std::cout<<"mask is empty!\n"; return {{}}; }
	bool allEmptyExcept5=true; //for if last plane is out of sync

	auto itmask = mask.begin();
//	std::cout<<"trigger "<<triggerNumberBegin<<" - "<<triggerNumberEnd<<"  at time "<<timestamp<<"\n";
	for (unsigned plane = 0; plane < mimosaHit->size(); plane++) {
//		std::cout<<"plane "<<plane<<" has "<<mimosaHit->at(plane).size()<<" hits, mask size "<<itmask->size()<<"\n";
		auto maskedHit = applyPixelMask(*itmask++, mimosaHit->at(plane));
		//convert hits to positions
		spaceHit.push_back(
				convertHits(maskedHit, plane, detector.planePosition[plane],
						detector.pixelsize, detector.pixelsize));
		if( (plane<=4 and not spaceHit.back().empty()) or (plane==5 and spaceHit.back().empty() ) ) allEmptyExcept5=false;
	}

	//Move entries one frame if last plane is out of sync
	if(not previousEntryHits.empty() and spaceHit.at(5).empty() ) {
		spaceHit[5]=previousEntryHits[5];
	}	else {
		if(allEmptyExcept5) {
			previousEntryHits=spaceHit;
			for(auto&s : spaceHit) s.clear();
		}
  	else previousEntryHits.clear();
	}


	//apply translation and rotation
	spaceHit=rotateAndShift(spaceHit);

	//set errors
	for(auto& hv: spaceHit) for(auto& h : hv) h.error=TVector3(detector.pixelsize/sqrt(12),detector.pixelsize/sqrt(12),1);

	return spaceHit;
}

std::vector<std::vector<PositionHit> >&  TelescopeTrackFitter::rotateAndShift(
		std::vector<std::vector<PositionHit> >& spaceHit) {
	//apply translation and rotation
	if (!shifts.empty())
		spaceHit = translateHits(spaceHit, shifts);
	if (!angles.empty())
		spaceHit = rotateHits(spaceHit, angles, hitsCentre);
	return spaceHit;
}

int TelescopeTrackFitter::getEntry(int iEvent) {
	if(iEvent>=nEvents) { std::cout<<"\nreached maximum number of events in telescope track fitter at event "<<iEvent<<"\n"; return false;}
	//get entry
	int nb=hitTable->GetEntry(iEvent);
	if(!nb) {
		std::cout<<"\nget entry returned 0 from telescope at event "<<iEvent<<"\n";
	}
	if (!mimosaHit) {
		cout << "Did not find mimosaHit!" << endl;
		throw 1;
		return false;
	}
	return nb;
}

void TelescopeTrackFitter::drawEvent(
		const std::vector<std::vector<PositionHit> >& spaceHit,
		const std::vector<FitResult3D>& fits) {
	HoughTransformer::drawClusters(spaceHit, detector);
	//			HoughTransformer::drawClusters(houghClusters, detector);

	for (auto& f : fits)
		f.draw( detector.zmin(), detector.zmax() );
	gPad->Update();
}

void TelescopeTrackFitter::fitTracks(std::string outputfilename) {
	residualHistograms=unique_ptr<ResidualHistogrammer>(new ResidualHistogrammer(outputfilename, detector));
	if(makeTrackHistograms) {
		trackHistograms=unique_ptr<TrackHistogrammer>(new TrackHistogrammer("") );
		trackHistogramsFirst=unique_ptr<TrackHistogrammer>(new TrackHistogrammer("First") );
		trackHistogramsSecond=unique_ptr<TrackHistogrammer>(new TrackHistogrammer("Second") );
	}
	if(generateHits) {
		simulatedQuadxResidual=unique_ptr<TH1D>(new TH1D("simulatedQuadxResidual","simulated quad x-residual;x-residual [mm];",100,-0.1,0.1));
		simulatedQuadyResidual=unique_ptr<TH1D>(new TH1D("simulatedQuadyResidual", "simulated quad y-residual; y-residual [mm];",100,-0.1,0.1));
	}

	//make fit three
	FitResult3D fit, fitFirstThree, fitSecondThree;
	TTree fitTree("fitTree", "fitTree");
	fitTree.Branch("fit", &fit);
	fitTree.Branch("fitFirstThree", &fitFirstThree);
	fitTree.Branch("fitSecondThree", &fitSecondThree);


	//for calculation of com, means and rotation
	std::vector<double> hitsXSum(detector.nPlanes), hitsYSum(detector.nPlanes);
	std::vector<double>	residualXSum(detector.nPlanes), residualYSum(detector.nPlanes);
	std::vector<double> rotationZSum(detector.nPlanes), rotationZWeightSum(detector.nPlanes);
	std::vector<int> nHits(detector.nPlanes);
	double slope1Sum=0, slope2Sum=0;

	//loop over all entries
	long int nPassed=0,nClusters=0;
	for( int iEvent=0; iEvent<1E5 //nEvents
	; iEvent++ ) {

		if(!(iEvent%10000))
			std::cout<<"event "<<iEvent<<"/"<<nEvents<<"\r"<<std::flush;

		//get entry
		if(!getEntry(iEvent) ) continue;

		//apply mask and convert
		TVector3 simulatedQuadPosition;
		std::vector<std::vector<PositionHit> > spaceHit = generateHits? generateSpaceHits(detector, simulatedQuadPosition) : getSpaceHits();

		std::vector<FitResult3D> fits;

		//check event
		if( !passEvent(spaceHit) ) continue;
		if(not nPassed) std::cout<<"first entry with hits at "<<iEvent<<"\n";
		++nPassed;

//		if(displayEvent) {
//			drawEvent(spaceHit, {} ) ;
//			if( processDrawSignals() ) break; //returns true if stopped
//		}

		//sum x and y here to sum all hits including noise

		//hough transform
		houghTransform.minClusterSize=5; //DEBUG! was 5
		auto houghClusters = doBinnedClustering ? binnedClustering(spaceHit) : houghTransform(spaceHit);

		if(!houghClusters.size()) {
			if(displayEvent) std::cout<<"no hough clusters found\n";
			continue;
		}
		//require at least nmin planes to be hit
		houghClusters.remove_if([this](const HoughTransformer::HitCluster& hc){return hc.getNPlanesHit()<nMinPlanesHit; });
		if(!houghClusters.size()) {
			if(displayEvent) std::cout<<"no hough clusters found with the minimum number of planes hit\n";
			continue;
		}

		//first fit on outer planes and translate inner hits
		for(auto& hitCluster : houghClusters) {

			//fit track
			if(hitCluster.size()<2 || hitCluster.recalculateNPlanesHit()<=1) {
				if(displayEvent) std::cout<<"not enough planes hit in hough cluster\n";
				continue;
			}

			hitCluster=clusterHitsPerPlane(hitCluster, detector);

			fit=regressionFit3d(hitCluster);
			if(!fit.isValid()) {
				cerr<<"fit not valid!"<<endl;
				continue;
			}

			//remove outliers
			hitCluster=calculateResiduals(hitCluster, fit);
			hitCluster=cutOnResiduals(hitCluster, maxResidual);

			//refit on planes with selection
			HoughTransformer::HitCluster selectedHits;
			std::copy_if(hitCluster.begin(), hitCluster.end(), std::back_inserter(selectedHits), selectHitForRefit );
			if((int) selectedHits.size()<nMinPlanesHit || selectedHits.recalculateNPlanesHit()<nMinPlanesHit) {
				if(displayEvent) std::cout<<"not enough planes hit in hough cluster after refit selection\n";
				continue;
			}

			if(constructLineParallelToZ) fit = makeLinesParallelToZ( selectedHits.front().position.x, selectedHits.front().position.y );
			else fit=regressionFit3d(selectedHits);
			if(!fit.isValid()) {cerr<<"fit not valid!"<<endl; continue;	}


			if(doPartialFits) {
				HoughTransformer::HitCluster hitsFirstThree;
				std::copy_if(selectedHits.begin(), selectedHits.end(), std::back_inserter(hitsFirstThree), [](const PositionHit& h){ return h.chip<=2;} );
				fitFirstThree=regressionFit3d(hitsFirstThree);

				HoughTransformer::HitCluster hitsSecondThree;
				std::copy_if(selectedHits.begin(), selectedHits.end(), std::back_inserter(hitsSecondThree), [](const PositionHit& h){ return h.chip>=3 and h.chip<=5;} );
				fitSecondThree=regressionFit3d(hitsSecondThree);
			}


			//sum x and y for calculation of centre of mass from track hits
			if(hitCluster.size()) for(auto& h : hitCluster) {
				hitsXSum[h.chip]+=h.position.x;
				hitsYSum[h.chip]+=h.position.y;
				++nHits[h.chip];
			}

			hitCluster=calculateResiduals(hitCluster, fit);

			residualHistograms->fill(hitCluster, hitsCentre);

			if(simulatedQuadxResidual) {
				simulatedQuadxResidual->Fill(fit.xAt(220)-simulatedQuadPosition.x());
				simulatedQuadyResidual->Fill(fit.yAt(220)-simulatedQuadPosition.y());
			}

//			SimpleDetectorConfiguration drawConfig(-0.05,0.05,-0.05,0.05,detector.zmin(), detector.zmax());
//			HoughTransformer::drawCluster(hitCluster, drawConfig);
//			if( processDrawSignals() ) return; //returns true if stopped

			//sum residuals
			for(auto& h : hitCluster) {
//				cout<<r.x<<" "<<r.y<<endl;
				residualXSum[h.chip]+=h.residual.x;
				residualYSum[h.chip]+=h.residual.y;

				auto& rotationPoint=hitsCentre[h.chip];
				double xc=rotationPoint.first, yc=rotationPoint.second; //x and y center
				double hx=h.position.x-xc, hy=h.position.y-yc;
				double phi=(hy*h.residual.x-hx*h.residual.y)/(hx*hx+hy*hy);
				double weight=hx*hx+hy*hy;
				rotationZSum[h.chip]+=phi*weight; rotationZWeightSum[h.chip]+=weight;
			}

			//sum fit slope
			slope1Sum+=fit.XZ.slope; slope2Sum+=fit.YZ.slope;

			fits.push_back(fit);
			if(makeTrackHistograms) {
				trackHistograms->fill(fit);
				trackHistogramsFirst->fill(fitFirstThree);
				trackHistogramsSecond->fill(fitSecondThree);

				fitTree.Fill();
			}

			//give fit errors
//			int i=0;
//			for(auto& planez : detector.planePosition ) {
//				cout<<"plane "<<i++<<": fit at "<<fit.xAt(planez)<<" +- "<<fit.XZ.errorAt(planez)<<", "<<fit.yAt(planez)<<" +- "<<fit.YZ.errorAt(planez)<<std::endl;
//			}

			++nClusters;
		}

		if(displayEvent) {
			drawEvent(spaceHit, fits) ;
			if( processDrawSignals() ) break; //returns true if stopped
		}

	}

	std::cout<<"\npassed: "<<nPassed<<"/"<<nEvents<<" with "<<nClusters<<" tracks \n\n";

	//Recalculate centre of mass of hits for each plane
	if(recalculateCOM) {
		cout<<"centre of chip is "<<detector.planexmax()/2<<", "<<detector.planeymax()/2<<endl;
		for(int i=0; i<detector.nPlanes; ++i) {
			hitsCentre[i]={ hitsXSum[i]/nHits[i], hitsYSum[i]/nHits[i] };
			cout<<"chip "<<i<<" hits centre of mass is: ("<<hitsCentre[i].first<<", "<<hitsCentre[i].second<<")"<<endl;
		}
		recalculateCOM=false;
	}

	for(int i=0; i<detector.nPlanes; i++) {
		cout<<"chip "<<i<<" average residual is "<<residualXSum[i]<<"/"<<nHits[i]<<", "<<residualYSum[i]<<"/"<<nHits[i]<<endl;
		averageResidualFromSum[i]= { residualXSum[i]/nHits[i], residualYSum[i]/nHits[i] };
		rotationZFromSum[i]= rotationZSum[i]/rotationZWeightSum[i];
	}

	slope1FromSum=slope1Sum/nClusters;
	slope2FromSum=slope2Sum/nClusters;
	cout<<"slopes are "<<slope1FromSum<<" "<<slope2FromSum<<endl;

	fitTree.Write();

//	cin.get();

}

std::vector<FitResult3D> TelescopeTrackFitter::getFits(std::vector<std::vector<PositionHit> >& spaceHit) {
	std::vector<FitResult3D> fits;
	partialFits.clear();

	//check event
	if( !passEvent(spaceHit) ) return fits;

	//sum x and y here to sum all hits including noise

	//hough transform
	auto houghClusters = houghTransform(spaceHit);

	//require at least nmin planes to be hit
//	const int nMinPlanesHit=4;
	houghClusters.remove_if([this](const HoughTransformer::HitCluster& hc){return hc.getNPlanesHit()<nMinPlanesHit; });
	if(!houghClusters.size()) {
		 return fits;
	}

	//first fit on outer planes and translate inner hits
	for(auto& hitCluster : houghClusters) {

		//fit track
		if(hitCluster.size()<2 || hitCluster.getNPlanesHit()<nMinPlanesHit) continue;
//			if(hitCluster.isPlaneHit(0)+hitCluster.isPlaneHit(1)+hitCluster.isPlaneHit(2) < 2) continue;
//			if(hitCluster.isPlaneHit(3)+hitCluster.isPlaneHit(4)+hitCluster.isPlaneHit(5) < 2) continue;

		auto fit=regressionFit3d(hitCluster);
		if(!fit.isValid()) {cerr<<"fit not valid!"<<endl;  continue;	}

		hitCluster=clusterHitsPerPlane(hitCluster, detector);

		//remove outliers
		hitCluster=calculateResiduals(hitCluster, fit);
		hitCluster=cutOnResiduals(hitCluster, maxResidual);

		//refit on planes with selection
		HoughTransformer::HitCluster selectedHits;
		std::copy_if(hitCluster.begin(), hitCluster.end(), std::back_inserter(selectedHits), selectHitForRefit );
		if(selectedHits.size()<2 || selectedHits.recalculateNPlanesHit()<nMinPlanesHit) continue;
//			if(selectedHits.isPlaneHit(0)+selectedHits.isPlaneHit(1)+selectedHits.isPlaneHit(2) < 2) continue;
//			if(selectedHits.isPlaneHit(3)+selectedHits.isPlaneHit(4)+selectedHits.isPlaneHit(5) < 2) continue;
		if(constructLineParallelToZ) fit = makeLinesParallelToZ( selectedHits.front().position.x, selectedHits.front().position.y );
		else fit=regressionFit3d(selectedHits);

		if(!fit.isValid()) {cerr<<"fit not valid!"<<endl; continue;	}

		//extra cut on scatter difference between first and last plane
		HoughTransformer::HitCluster hitsFirstThree, hitsSecondThree;
		std::copy_if(selectedHits.begin(), selectedHits.end(), std::back_inserter(hitsFirstThree), [](const PositionHit& h){ return h.chip<=2;} );
		std::copy_if(selectedHits.begin(), selectedHits.end(), std::back_inserter(hitsSecondThree), [](const PositionHit& h){ return h.chip>=3; } );
		if(hitsFirstThree.recalculateNPlanesHit()<=1 or hitsSecondThree.recalculateNPlanesHit()<=1) continue;
		auto fitFirst=regressionFit3d(hitsFirstThree);
		auto fitSecond=regressionFit3d(hitsSecondThree);
		if( fabs(fitFirst.XZ.slope-fitSecond.XZ.slope) > 1E-3 or fabs(fitFirst.YZ.slope-fitSecond.YZ.slope) > 1E-3 ) continue;

		hitCluster=calculateResiduals(hitCluster, fit);
		fits.push_back(fit);
		partialFits.push_back( {fitFirst,fitSecond} );
	}

	return fits;
}

std::vector<std::pair<double, double> > TelescopeTrackFitter::getMeanResiduals() {
//	auto means= residualHistograms->getMeansOfPlanes();
	auto means= averageResidualFromSum;
	for(auto& m : means ) cout<<"shift ("<<m.first<<", "<<m.second<<")"<<endl;
	return means;
}

std::vector<double> TelescopeTrackFitter::getRotations() {
//	auto rotation= residualHistograms->getRotationOfPlanes();
	auto rotation= rotationZFromSum;
	int i=0;
	for(auto& r : rotation ) cout<<"chip "<<i++<<" rotation: "<<r<<" = "<<r/M_PI*180.<<endl;
	return rotation;
}

void TelescopeTrackFitter::setShifts(
		const std::vector<std::pair<double, double> >& shiftsIn) {
		shifts=shiftsIn;
}

void TelescopeTrackFitter::setAngles(const std::vector<double>& anglesIn) {
	angles=anglesIn;
}

void TelescopeTrackFitter::addToShifts(
		const std::vector<std::pair<double, double> >& shiftsExtra) {
	if(shifts.empty()) return setShifts(shiftsExtra);

	auto itExtra= shiftsExtra.begin();
	for(auto& s : shifts) {
		s.first+=itExtra->first;
		s.second+=itExtra->second;
		++itExtra;
	}
}

void TelescopeTrackFitter::addToAngles(const std::vector<double>& anglesExtra) {
	if(angles.empty()) return setAngles(anglesExtra);

	auto itExtra=anglesExtra.begin();
	for(auto& a : angles) {
		a+=*itExtra++;
	}
}

const std::vector<std::pair<double, double> >& TelescopeTrackFitter::getShifts() const {
	return shifts;
}

const std::vector<double>& TelescopeTrackFitter::getAngles() const {
	return angles;
}

std::pair<double, double> TelescopeTrackFitter::getSlopes() const {
	return {slope1FromSum, slope2FromSum};
}

void TelescopeTrackFitter::setSlopes(std::pair<double, double> slopes) {
	houghTransform.angleOfTracksX=slopes.first;
	houghTransform.angleOfTracksY=slopes.second;
}

void TelescopeTrackFitter::saveAlignment(std::string file) {

	//load alignment parameter, because we need to save the other parameters in the same file
	Alignment alignment(file);

	alignment.mimosa.shifts=shifts;
	alignment.mimosa.COMs=hitsCentre;
	alignment.mimosa.angles=angles;
	alignment.mimosa.slopes={houghTransform.angleOfTracksX, houghTransform.angleOfTracksY};

	alignment.write(file);
}


void TelescopeTrackFitter::setCentres(
		const std::vector<std::pair<double, double> >& COMs) {
	hitsCentre=COMs;
}

void TelescopeTrackFitter::setAlignment(const Alignment& alignment) {
	setShifts( alignment.mimosa.shifts ); //xy shift of planes
	setAngles( alignment.mimosa.angles ); //planes rotation along rz-axis in xy plane
	setCentres( alignment.mimosa.COMs );
	setSlopes( alignment.mimosa.slopes );
}
