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

#include "../Alignment/Alignment.h"
#include "DetectorConfiguration.h"
#include "/user/cligtenb/rootmacros/getObjectFromFile.h"
//#include "../../rootmacros/getObjectFromFile.h"
#include "ResidualHistogrammer.h"
#include "transformHits.h"
#include "makeNoisyPixelMask.cpp"

using namespace std;

TelescopeTrackFitter::TelescopeTrackFitter(std::string inputfile, const TelescopeConfiguration& detector) :
	detector(detector),
	houghTransform(detector.xmin(), detector.xmax(), detector.ymin(), detector.ymax(), 200/*xbins*/, 500 /*ybins*/ ),
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
	nEvents=hitTable->GetEntriesFast();
	//    unsigned short triggerNumberBegin, triggerNumberEnd;
	hitTable->SetBranchAddress("triggerNumberBegin", &triggerNumberBegin);
	hitTable->SetBranchAddress("triggerNumberEnd", &triggerNumberEnd);
	hitTable->SetBranchAddress("timestamp", &timestamp);

}

TelescopeTrackFitter::~TelescopeTrackFitter() {
	file->Close();
	residualHistograms=nullptr; //make sure to close residualHistograms before track, because they use the same file
	trackHistograms=nullptr; //todo: combine both in a proper way;
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
	if(iEvent>=nEvents) return false;
	//get entry
	int nb=hitTable->GetEntry(iEvent);
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
	if(makeTrackHistograms) trackHistograms=unique_ptr<TrackHistogrammer>(new TrackHistogrammer(detector) );

	//for calculation of com, means and rotation
	std::vector<double> hitsXSum(detector.nPlanes), hitsYSum(detector.nPlanes);
	std::vector<double>	residualXSum(detector.nPlanes), residualYSum(detector.nPlanes);
	std::vector<double> rotationZSum(detector.nPlanes), rotationZWeightSum(detector.nPlanes);
	std::vector<int> nHits(detector.nPlanes);
	double slope1Sum=0, slope2Sum=0;

	//loop over all entries
	long int nPassed=0,nClusters=0;
	for( int iEvent=0; iEvent<nEvents; iEvent++ ) {

		if(!(iEvent%10000))
			std::cout<<"event "<<iEvent<<"/"<<nEvents<<"\r"<<std::flush;

		//get entry
		if(!getEntry(iEvent) ) continue;

		//apply mask and convert
		std::vector<std::vector<PositionHit> > spaceHit = getSpaceHits();

		std::vector<FitResult3D> fits;

		//check event
		if( !passEvent(spaceHit) ) continue;
		if(not nPassed) std::cout<<"first entry with hits at "<<iEvent<<"\n";
		++nPassed;

		if(displayEvent) {
			drawEvent(spaceHit, {} ) ;
			if( processDrawSignals() ) break; //returns true if stopped
		}

		//sum x and y here to sum all hits including noise

		//hough transform
		auto houghClusters = doBinnedClustering ? binnedClustering(spaceHit) : houghTransform(spaceHit);

		//require at least nmin planes to be hit
		const int nMinPlanesHit=4;
		houghClusters.remove_if([](const HoughTransformer::HitCluster& hc){return hc.getNPlanesHit()<nMinPlanesHit; });
		if(!houghClusters.size()) {
			continue;
		}

		//first fit on outer planes and translate inner hits
		for(auto& hitCluster : houghClusters) {

			//fit track
			if(hitCluster.size()<2 || hitCluster.getNPlanesHit()<=1) continue;
//			if(hitCluster.isPlaneHit(0)+hitCluster.isPlaneHit(1)+hitCluster.isPlaneHit(2) < 2) continue;
//			if(hitCluster.isPlaneHit(3)+hitCluster.isPlaneHit(4)+hitCluster.isPlaneHit(5) < 2) continue;

			auto fit=regressionFit3d(hitCluster);
			if(!fit.isValid()) {cerr<<"fit not valid!"<<endl;  continue;	}

			//remove outliers
			hitCluster=calculateResiduals(hitCluster, fit);
			hitCluster=cutOnResiduals(hitCluster, maxResidual);

			//refit on planes with selection
			HoughTransformer::HitCluster selectedHits;
			std::copy_if(hitCluster.begin(), hitCluster.end(), std::back_inserter(selectedHits), selectHitForRefit );
			if(selectedHits.size()<2 || selectedHits.recalculateNPlanesHit()<=1) continue;
//			if(selectedHits.isPlaneHit(0)+selectedHits.isPlaneHit(1)+selectedHits.isPlaneHit(2) < 2) continue;
//			if(selectedHits.isPlaneHit(3)+selectedHits.isPlaneHit(4)+selectedHits.isPlaneHit(5) < 2) continue;
			if(constructLineParallelToZ) fit = makeLinesParallelToZ( selectedHits.front().position.x, selectedHits.front().position.y );
			else fit=regressionFit3d(selectedHits);

			if(!fit.isValid()) {cerr<<"fit not valid!"<<endl; continue;	}

			//sum x and y for calculation of centre of mass from track hits
			if(hitCluster.size()) for(auto& h : hitCluster) {
				hitsXSum[h.chip]+=h.position.x;
				hitsYSum[h.chip]+=h.position.y;
				++nHits[h.chip];
			}

			hitCluster=calculateResiduals(hitCluster, fit);

			residualHistograms->fill(hitCluster, hitsCentre);

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
			if(makeTrackHistograms) { trackHistograms->fill(fit); }

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

	std::cout<<"\npassed: "<<nPassed<<"/"<<nEvents<<" with "<<nClusters<<"tracks \n";

	//Recalculate centre of mass of hits for each plane
	if(recalculateCOM) {
		cout<<"centre of chip is "<<detector.planexmax()/2<<", "<<detector.planeymax()/2<<endl;
		for(int i=0; i<detector.nPlanes; ++i) {
			hitsCentre[i]={ hitsXSum[i]/nHits[i], hitsYSum[i]/nHits[i] };
			cout<<"Hits centre of mass is: ("<<hitsCentre[i].first<<", "<<hitsCentre[i].second<<")"<<endl;
		}
		recalculateCOM=false;
	}

	for(int i=0; i<detector.nPlanes; i++) {
		cout<<"average residual is "<<residualXSum[i]<<"/"<<nHits[i]<<", "<<residualYSum[i]<<"/"<<nHits[i]<<endl;
		averageResidualFromSum[i]= { residualXSum[i]/nHits[i], residualYSum[i]/nHits[i] };
		rotationZFromSum[i]= rotationZSum[i]/rotationZWeightSum[i];
	}

	slope1FromSum=slope1Sum/nClusters;
	slope2FromSum=slope2Sum/nClusters;
	cout<<"slopes are "<<slope1FromSum<<" "<<slope2FromSum<<endl;

//	cin.get();

}

std::vector<FitResult3D> TelescopeTrackFitter::getFits(std::vector<std::vector<PositionHit> >& spaceHit) {
	std::vector<FitResult3D> fits;

	//check event
	if( !passEvent(spaceHit) ) return fits;

	//sum x and y here to sum all hits including noise

	//hough transform
	auto houghClusters = houghTransform(spaceHit);

	//require at least nmin planes to be hit
	const int nMinPlanesHit=4;
	houghClusters.remove_if([](const HoughTransformer::HitCluster& hc){return hc.getNPlanesHit()<nMinPlanesHit; });
	if(!houghClusters.size()) {
		 return fits;
	}

	//first fit on outer planes and translate inner hits
	for(auto& hitCluster : houghClusters) {

		//fit track
		if(hitCluster.size()<2 || hitCluster.getNPlanesHit()<=1) continue;
//			if(hitCluster.isPlaneHit(0)+hitCluster.isPlaneHit(1)+hitCluster.isPlaneHit(2) < 2) continue;
//			if(hitCluster.isPlaneHit(3)+hitCluster.isPlaneHit(4)+hitCluster.isPlaneHit(5) < 2) continue;

		auto fit=regressionFit3d(hitCluster);
		if(!fit.isValid()) {cerr<<"fit not valid!"<<endl;  continue;	}

		//remove outliers
		hitCluster=calculateResiduals(hitCluster, fit);
		hitCluster=cutOnResiduals(hitCluster, maxResidual);

		//refit on planes with selection
		HoughTransformer::HitCluster selectedHits;
		std::copy_if(hitCluster.begin(), hitCluster.end(), std::back_inserter(selectedHits), selectHitForRefit );
		if(selectedHits.size()<2 || selectedHits.recalculateNPlanesHit()<=1) continue;
//			if(selectedHits.isPlaneHit(0)+selectedHits.isPlaneHit(1)+selectedHits.isPlaneHit(2) < 2) continue;
//			if(selectedHits.isPlaneHit(3)+selectedHits.isPlaneHit(4)+selectedHits.isPlaneHit(5) < 2) continue;
		if(constructLineParallelToZ) fit = makeLinesParallelToZ( selectedHits.front().position.x, selectedHits.front().position.y );
		else fit=regressionFit3d(selectedHits);

		if(!fit.isValid()) {cerr<<"fit not valid!"<<endl; continue;	}


		hitCluster=calculateResiduals(hitCluster, fit);
		fits.push_back(fit);
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
	for(auto& r : rotation ) cout<<"rotation: "<<r<<" = "<<r/M_PI*180.<<endl;
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
