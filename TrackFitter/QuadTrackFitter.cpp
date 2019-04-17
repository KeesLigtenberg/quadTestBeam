/*
 * LaserDataFitter.cpp
 *
 *  Created on: Jul 16, 2018
 *      Author: cligtenb
 */

#include "QuadTrackFitter.h"

#include <algorithm>
#include <array>

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TTree.h"
#include "TPolyLine3D.h"

#include "PositionHit.h"
#include "../Alignment/TimeWalkCorrector.h"
#include "../TelescopeTrackFitter/HoughTransformer.h"

#pragma link C++ class std::vector<int>+;

namespace {

	const int nChips=4;
//	struct shift {double x,y;};
//	std::array<shift,nChips> chipShifts={{ {24.1,15.7}, {38.0, 26.3}, {23.7,26.3}, {9.9,15.7} }};


	std::vector<PositionHit> convertHitsQuad( std::vector<Hit>* chips [], double driftSpeed, double t0Offset=0 /*value to subtract from time*/) {
		std::vector<PositionHit> hits;
		for(int i=0; i<nChips; i++) {
			auto chipHits=convertHitsTPC(*chips[i], i, driftSpeed, t0Offset); //t0 is subtracted from hit time
			hits.insert(hits.end(), chipHits.begin(), chipHits.end() );
		}
		return hits;
	}

	//selection functions for laser points
	bool selectByRectangularCut(const Vec3& v) {
		struct {	double x,y;	} chipCenters[4] = { {17,8.5}, {31,8.5}, {17,33.5}, {31,33.5} } ;
		double maxDistance=5.1;
		for(auto& center : chipCenters) {
			if(std::fabs(center.x-v.x)<maxDistance and std::fabs(center.y-v.y)<maxDistance) {
				return true;
			}
		}
		return false;
	}
	struct cutByMinimumDistanceFromEdge{
		const Alignment& alignment;
		double minDistance;//mm
		bool operator() (const Vec3& laser) {
			if(std::any_of(alignment.chips.begin(), alignment.chips.end(), [&laser, this](const ChipAlignment& ca){
					return ca.getDistanceFromEdge(laser)<minDistance;
				})) {
				return false;
			}
			return true;
		}
	};


	void drawChipFromCorners(std::array<TVector3, 4> corners) {
		const int npoints=5;
		double x[npoints] = {};
		double y[npoints] = {};
		double z[npoints] = {};
		for(int i=0; i<npoints; i++) {
			x[i]=corners[i%4].x();
			y[i]=corners[i%4].y();
			z[i]=corners[i%4].z();
		}
		TPolyLine3D l( npoints, x, y, z );
		l.SetLineColor(kAzure-8);
		l.SetLineWidth(2);
		l.DrawClone();
	}

	void drawQuadOutline(const Alignment& align, double z) {
		for(int i=0; i<nChips; i++) {
			auto corners=align.getChipCorners(i);
			for(auto& c: corners) c.SetZ( z );
			drawChipFromCorners(corners);
		}
	}


	struct PositionRange { double xmin, xmax, ymin, ymax; };
	inline PositionRange rangeFromCorners(const std::vector<TVector3>& corners) {
		auto xrange=std::minmax_element(corners.begin(), corners.end(), [](const TVector3& v, const TVector3& w){return v.x()<w.x();});
		auto yrange=std::minmax_element(corners.begin(), corners.end(), [](const TVector3& v, const TVector3& w){return v.y()<w.y();});
		auto xmin=xrange.first->x(),xmax=xrange.second->x(),ymin=yrange.first->y(),ymax=yrange.second->y();
		return {xmin,xmax,ymin,ymax};
	}


	PositionHit& correctzPerPixel(PositionHit& h){
//		static TProfile2D* zResidualPerPixel=getObjectFromFile<TProfile2D>("zResidualByPixel", "pixelCorrection.root");
		static TProfile2D* zResidualPerPixel=getObjectFromFile<TProfile2D>("correction", "/project/lepcol/users/cligtenb/quadTestBeam/run668/pixelCorrectionSmaller.root");
		static TProfile* rowCorrection=getObjectFromFile<TProfile>("rowzCorrection", "/project/lepcol/users/cligtenb/quadTestBeam/run668/pixelCorrectionSmaller.root");
		static TProfile* columnCorrection=getObjectFromFile<TProfile>("columnzCorrection", "/project/lepcol/users/cligtenb/quadTestBeam/run668/pixelCorrectionSmaller.root");
		auto bin=zResidualPerPixel->GetBin((h.column)%zResidualPerPixel->GetNbinsX()+1,(h.row)%zResidualPerPixel->GetNbinsY()+1);
//		if(zResidualPerPixel->GetBinEntries(bin) > 50 )
				h.position.z+=zResidualPerPixel->GetBinContent(bin);
				h.position.z+=std::max(std::min(rowCorrection->GetBinContent(h.row+1), 0.1),-0.1) ;
				h.position.z+=std::max(std::min(columnCorrection->GetBinContent(h.column+1), 0.1),-0.1) ;
		return h;
	}

}




struct ChipHistogrammer {
	ChipHistogrammer(const char* name, const Alignment& align) : ChipHistogrammer(std::string(name), align) {}
	ChipHistogrammer(std::string name, const Alignment& align) ;

	std::unique_ptr<TH1I> nHits;
	std::unique_ptr<TH1D> driftTime, ToT, nShiftedTrigger;
	std::unique_ptr<TH2D> pixelHitMap;
	std::unique_ptr<TH1D> xRotation, yRotation, zRotation;

	std::unique_ptr<TH2D> zResidualByToT, xResidualByToT, zResidualByToTCorrected,zResidualByDriftTime;

	struct residualHistograms {
		residualHistograms(std::string name) : name{name} {};
		std::string name;
		std::unique_ptr<TH1D> xResidual{}, yResidual{}, zResidual{};
		std::unique_ptr<TH1D> zHit{};
		std::unique_ptr<TH2D> xResidualByz{}, yResidualByz{},zResidualByz{};
		std::unique_ptr<TH3D> zResidualByzByToT{};
		std::unique_ptr<TH2D>	positionHitMap{};
		std::unique_ptr<TProfile2D> xResidualByPosition{}, zResidualByPosition{};
		std::unique_ptr<TProfile2D> xResidualByXZ{}, zResidualByXZ{};
		std::vector< std::unique_ptr<TProfile2D> > xResidualByXYZ{};
		void fill(const Vec3& position, const Vec3& residual, double ToT) {
			xResidual->Fill(residual.x);
			yResidual->Fill(residual.y);
			zResidual->Fill(residual.z);
			xResidualByz->Fill(position.z,residual.x);
			yResidualByz->Fill(position.z,residual.y);
			zResidualByz->Fill(position.z,residual.z);
			zResidualByzByToT->Fill(position.z,residual.z, ToT*.025);
			positionHitMap->Fill(position.x, position.y);
			zHit->Fill(position.z);
			xResidualByPosition->Fill(position.x,position.y,residual.x);
			zResidualByPosition->Fill(position.x,position.y,residual.z);

			int zbin=(position.z+1)/2;
			if(zbin>=0 and zbin<=4) xResidualByXYZ[zbin]->Fill(position.x,position.y, residual.x);
			xResidualByXZ->Fill(position.x,position.z,residual.x);
			zResidualByXZ->Fill(position.x,position.z,residual.z);
		}
		void createHistograms(const PositionRange& range);
	} global{"_global"}, local{"_local"}, locExp{"_locExp"}, corrected{"_corrected"}, tighterCuts{"_tighterCuts"}, sideBands{"_sideBands"};
	std::unique_ptr<TProfile2D> xResidualByPixel, zResidualByPixel;

	void fillHit(const PositionHit& h);
	void fillRotation(const TVector3& pos, const TVector3& residual, const TVector3& COM);
	void fillRotation(const PositionHit& h, const TVector3& COM) { fillRotation(h.position, h.residual, COM); }
	void fillEvent();
	void fillTimewalk(double localZResidual, double ToT, const TimeWalkCorrector& twc) {
		zResidualByToT->Fill(ToT*25E-3, localZResidual+twc.getCorrection(ToT));
	}
private:
	std::string dirName;
	int hitsCounter=0;
};

ChipHistogrammer::ChipHistogrammer(std::string name, const Alignment& align) : dirName(name) {
	using namespace std;

	//change directory
	auto startDir=gDirectory;
	gDirectory->mkdir( dirName.c_str() )->cd();

	nHits=		unique_ptr<TH1I>(new TH1I{"nHits", "Number of hits per event;Number of hits;Entries", 200,0,400});
	constexpr double ToABinWidth=1.5625E-3;
	driftTime=unique_ptr<TH1D>(new TH1D{"driftTime","Drift time;Drift time [#mus];Hits", int(0.8/ToABinWidth),-0.39999,0.4}); //zero would be rounded up (wrong), so move bin slightly to 0.099999
	ToT=			unique_ptr<TH1D>(new TH1D{"ToT", "Time over threshold;ToT [#mus];Hits",80,0,2});
	nShiftedTrigger= unique_ptr<TH1D>(new TH1D{"nShiftedTrigger", "nShiftedTrigger;Shifted time [409.6 #mus]; Entries", 500,0,500});

	pixelHitMap=unique_ptr<TH2D>(new TH2D{"pixelHitMap", "Hitmap by pixel;Columns;Rows", 256,0, 256, 256,0,256});

	xRotation=unique_ptr<TH1D>(new TH1D{"xRotation", "x-rotation;x-rotation [rad.]; Hits", 100,-0.5,0.5});
	yRotation=unique_ptr<TH1D>(new TH1D{"yRotation", "y-rotation;y-rotation [rad.]; Hits", 100,-0.5,0.5});
	zRotation=unique_ptr<TH1D>(new TH1D{"zRotation", "z-rotation;z-rotation [rad.]; Hits", 100,-0.5,0.5});


	zResidualByToT=std::unique_ptr<TH2D>(new TH2D{"zResidualByToT", "z-residual by ToT;ToT [#mus]; z-residual [mm]", 100,0,2.5, 200,-5,5});
	xResidualByToT=std::unique_ptr<TH2D>(new TH2D{"xResidualByToT", "x-residual by ToT;ToT [#mus]; z-residual [mm]", 100,0,2.5, 200,-5,5});
	zResidualByToTCorrected=std::unique_ptr<TH2D>(new TH2D{"zResidualByToTCorrected", "z-residual by ToT;ToT [#mus]; z-residual [mm]", 100,0,2.5, 200,-5,5});
	zResidualByDriftTime=std::unique_ptr<TH2D>(new TH2D{"zResidualByDriftTime", "z-residuals as a function of drift time;Drift time [#mus];z-residual [mm]", int(0.8/ToABinWidth),-0.39999,0.4,50,-2,2});

	auto rangeGlobal=rangeFromCorners( align.getAllChipCorners() );
//	auto rangeQuad=rangeFromCorners( align.getAllChipCornersQuad() ); //automatic range
	PositionRange rangeQuad{0,29.04,-14.55,25.05};

	global.createHistograms(rangeGlobal);
	local.createHistograms(rangeQuad);
	locExp.createHistograms(rangeQuad);
	corrected.createHistograms(rangeQuad);
	tighterCuts.createHistograms(rangeQuad);
	sideBands.createHistograms(rangeQuad);

	xResidualByPixel=std::unique_ptr<TProfile2D>(new TProfile2D{"xResidualByPixel", "mean x-residual by pixel;Columns;Rows", 256,0,256,256,0,256});
	zResidualByPixel=std::unique_ptr<TProfile2D>(new TProfile2D{"zResidualByPixel", "mean z-residual by pixel;Columns;Rows", 256,0,256,256,0,256});

	xResidualByPixel->SetMinimum(-0.1), xResidualByPixel->SetMaximum(0.1);
	zResidualByPixel->SetMinimum(-0.1), zResidualByPixel->SetMaximum(0.1);

	startDir->cd();
}

void ChipHistogrammer::residualHistograms::createHistograms(const PositionRange& range) {
	using namespace std;
	const int nBins=160, nBinsZ=150;

	auto startDir=gDirectory;
	if(!startDir) std::cout<<"current directory is not valid!\n";
	auto newDir=startDir->mkdir( name.substr(1,std::string::npos).c_str() );
	if(!newDir) std::cout<<"could not make directory: " <<name<<"\n";
	newDir->cd();

	xResidual=unique_ptr<TH1D>(new TH1D{("xResidual"), "x-residual;x-residual [mm]; Hits", nBins,-2,2});
	yResidual=unique_ptr<TH1D>(new TH1D{("yResidual"), "y-residual;y-residual [mm]; Hits", nBins,-2,2});
	zResidual=unique_ptr<TH1D>(new TH1D{("zResidual"), "z-residual;z-residual [mm]; Hits", nBins,-4,4});

	xResidualByz=std::unique_ptr<TH2D>(new TH2D{("xResidualByz"), "x-residuals as a function of drift distance;Drift distance [mm];x-residual [mm]", nBinsZ,-1,14,nBins,-2,2});
	yResidualByz=std::unique_ptr<TH2D>(new TH2D{("yResidualByz"), "y-residuals as a function of drift distance;Drift distance [mm];y-residual [mm]", nBinsZ,-1,14,nBins,-2,2});
	zResidualByz=std::unique_ptr<TH2D>(new TH2D{("zResidualByz"), "z-residuals as a function of drift distance;Drift distance [mm];z-residual [mm]", nBinsZ,-1,14,nBins,-2,2});
	zResidualByzByToT=std::unique_ptr<TH3D>(new TH3D{("zResidualByzByToT"), "z-residuals as a function of drift distance;Drift distance [mm];z-residual [mm];ToT [#mus]", nBinsZ,-1,14,nBins,-2,2,100,0,2.5});

	double zmin=-5, zmax=15;
	zHit = unique_ptr<TH1D>(new TH1D{ ("zHit"), "z-position of hits;z [mm]; Hits", 200,zmin,zmax });
	double xmax=range.xmax, xmin=range.xmin, ymin=range.ymin, ymax=range.ymax;
	positionHitMap=unique_ptr<TH2D>(new TH2D{ ("positionHitMap"), "Hitmap with positions;x-position [mm];y-position [mm]",
		int((xmax-xmin)/.055),xmin,xmax,int((ymax-ymin)/.055),ymin, ymax});

	xResidualByPosition=std::unique_ptr<TProfile2D>(new TProfile2D{ ("xResidualByPosition"), "mean x-residual by position;x-position [mm]; y-position [mm]; mean x-residual [mm]",
		int((xmax-xmin)/.055),xmin,xmax,int((ymax-ymin)/.055),ymin, ymax});
	zResidualByPosition=std::unique_ptr<TProfile2D>(new TProfile2D{ ("zResidualByPosition"), "mean z-residual by position;x-position [mm]; y-position [mm]; mean z-residual [mm]",
		int((xmax-xmin)/.055),xmin,xmax,int((ymax-ymin)/.055),ymin, ymax});
	for(auto prof2d : {&xResidualByPosition, &zResidualByPosition} ) {
		(*prof2d)->SetMinimum(-0.1), (*prof2d)->SetMaximum(0.1);
	}

	for(int i=0; i<5; i++) {
		xResidualByXYZ.emplace_back( new TProfile2D{ ("xResidualByXYZ"+to_string(i)).c_str(),
			("From z "+std::to_string(i*2-1)+" to "+std::to_string((i+1)*2-1)+";x [mm]; y[mm];z [mm]; mean x-residual [mm]").c_str(),
			int((xmax-xmin)/.055),xmin,xmax,int((ymax-ymin)/.055),ymin, ymax} );
	}

	xResidualByXZ=std::unique_ptr<TProfile2D>(new TProfile2D{ ("xResidualByXZ"), "mean x-residual by XZ;x-position [mm]; y-position [mm]; mean x-residual [mm]",
		int((xmax-xmin)/.055),xmin,xmax,100,zmin, zmax});
	zResidualByXZ=std::unique_ptr<TProfile2D>(new TProfile2D{ ("zResidualByXZ"), "mean z-residual by XZ;x-position [mm]; y-position [mm]; mean z-residual [mm]",
		int((xmax-xmin)/.055),xmin,xmax,100,zmin, zmax});
	for(auto prof2d : {&xResidualByXZ, &zResidualByXZ} ) {
		(*prof2d)->SetMinimum(-0.2), (*prof2d)->SetMaximum(0.2);
	}

	startDir->cd();
}

void ChipHistogrammer::fillHit(const PositionHit& h) {
	driftTime->Fill(h.driftTime/4096.*25E-3);
	ToT->Fill(h.ToT*25E-3);
	pixelHitMap->Fill(h.column, h.row);
	zResidualByToTCorrected->Fill(h.ToT*25E-3, h.residual.z);
	xResidualByToT->Fill(h.ToT*25E-3, h.residual.x);
	zResidualByDriftTime->Fill(h.driftTime/4096.*25E-3,h.residual.z);
	nShiftedTrigger->Fill(h.nShiftedTrigger);
	global.fill(h.position, h.residual, h.ToT);
	xResidualByPixel->Fill(h.column,h.row, h.residual.x);
	zResidualByPixel->Fill(h.column,h.row, h.residual.z);

	hitsCounter++;
}

void ChipHistogrammer::fillRotation(const TVector3& position, const TVector3& residual, const TVector3& COM) { //h in timepix coordinates, COM
	//todo: do proper 3d rotation alignment
	auto d=position-COM;
	//x
	double perp2x=d.Perp2({1,0,0});
	double phix=(d.z()*residual.y()-d.y()*residual.z())/perp2x;
	xRotation->Fill(phix, perp2x);

	//y
	double perp2y=d.Perp2({0,1,0});
	double phiy=(d.x()*residual.z()-d.z()*residual.x())/perp2y;
	yRotation->Fill(phiy, perp2y);
	//for the moment, until error are included, use only x residuals in y rotation!
//	double phiy=-h.residual.x/d.z();
//	yRotation->Fill(phiy,d.z()*d.z());

	//z
//	double phi=(d.y()*h.residual.x-d.x()*h.residual.y)/d.Perp2();
//	double weight=d.Perp2();
	double phi=residual.x()/d.y();
	double weight=d.y()*d.y();
	zRotation->Fill(phi,weight);
}
void ChipHistogrammer::fillEvent() {
	if(!hitsCounter) return;
	nHits->Fill(hitsCounter);
	hitsCounter=0;
}


QuadTrackFitter::QuadTrackFitter(std::string fileName) :
		reader("data", fileName) {}

QuadTrackFitter::~QuadTrackFitter() {
}

std::vector<PositionHit> QuadTrackFitter::getSpaceHits(const Alignment& alignment) {
	//retrieve hits
	const auto driftSpeed = alignment.driftSpeed.value / 4096 * 25; //in units of mm/ticks
	const auto t0Offset = alignment.t0Offset.value /25. * 4096; //convert from ns to ticks
	posHits = convertHitsQuad(reader.chip, driftSpeed, t0Offset);
	for (auto& h : posHits) {
		//apply alignment
		h = alignment.totCorrections[h.chip].correct(h);

		//corrections:

		//do bin correction for columns 192, 193 (new quad)
		if(h.column==192 or h.column==193) {
			const double toaCorrection=25*alignment.driftSpeed.value;
			h.position.z-=toaCorrection;
			h.driftTime-=toaCorrection/alignment.driftSpeed.value;
		}
		//do odd/even bin column correction
		const double oddEvenCorrection=0.02;
		if(h.column%2) {//odd
			h.position.z-=oddEvenCorrection;
		} else { //even
			h.position.z+=oddEvenCorrection;
		}
		//per pixel correction
		h=correctzPerPixel(h);

		h = alignment.timeWalks[h.chip].correct(h);
		h = alignment.transform(h);
		h.error = alignment.hitErrors.hitError( alignment.quad.rotateAndShiftBack(h.position).z );
	}
	return posHits;
}

void QuadTrackFitter::Loop(std::string outputFile,const Alignment& alignment) {

	TFile file(outputFile.c_str(), "RECREATE");

	Vec3 averageHitPosition;
	std::vector<int> nHits(4);
	std::vector<int> nHitsPassed(4);
	int nHitsPassedTotal{0};

	TTree fitResults("fitResults", "fitResults");
	fitResults.Branch("hits", "std::vector<PositionHit>", &posHits);
	fitResults.Branch("hitAverage", "Vec3", &averageHitPosition);
	fitResults.Branch("nHitsChip", "std::vector<int>", &nHits);
	fitResults.Branch("nHitsPassedChip", "std::vector<int>", &nHitsPassed);
	fitResults.Branch("nHitsPassed", &nHitsPassedTotal);

  auto zResidualByToT=std::unique_ptr<TH2D>(new TH2D{"zResidualByToT", "z-residual by ToT;ToT [#mus]; z-residual [mm]", 100,0,2.5, 200,-4,6});
  auto zResidualByToTCorrected=std::unique_ptr<TH2D>(new TH2D{"zResidualByToTCorrected", "z-residual by ToT;ToT [#mus]; z-residual [mm]", 100,0,2.5, 200,-4,6});

	std::array<ChipHistogrammer, nChips> hists{{
		{"chip1", alignment},
		{"chip2", alignment},
		{"chip3", alignment},
		{"chip4", alignment} }};
	ChipHistogrammer quad{"quad", alignment};

	auto nEntries=reader.tree->GetEntries();
	std::cout<<nEntries<<" entries\n";
	long long cEntry=0;
	while( getEntry(cEntry++) ) {
		if( !(cEntry%10000) ) std::cout<<"entry "<<cEntry<<" / "<<nEntries<<"\n";

		if(cEntry>=1E5) break;

		//retrieve hits
		posHits=getSpaceHits(alignment);
		nHits=countHitsPerChip(posHits,false);
		nHitsPassed=countHitsPerChip(posHits,true);
		nHitsPassedTotal=std::accumulate(nHitsPassed.begin(), nHitsPassed.end(),0);

//		if( (nHits[0] + nHits[1]>70) ) {
//			for(auto& h : posHits) {
//				if(h.chip==2 or h.chip==3 )
//					h.flag=PositionHit::Flag::debug;
//			}
//		} else if ( (nHits[2]+nHits[3]>70) ) {
//			for(auto& h : posHits) {
//				if(h.chip==0 or h.chip==1 )
//					h.flag=PositionHit::Flag::debug;
//			}
//		} else {
//			continue;
//		}

		averageHitPosition={0,0,0};
		for(auto& h : posHits) {

			if(h.flag==PositionHit::Flag::highResidualxy) continue;
			hists[h.chip].fillTimewalk(h.ToT, h.residual.z, alignment.timeWalks[h.chip]);
			quad.fillTimewalk(h.ToT, h.residual.z, alignment.timeWalks[h.chip]);
			zResidualByToTCorrected->Fill(h.ToT*25E-3,h.residual.z);
			zResidualByToT->Fill(h.ToT*25E-3,h.residual.z+alignment.timeWalks[h.chip].getCorrection(h.ToT));

			if(h.flag!=PositionHit::Flag::valid) continue;

			hists[h.chip].fillHit(h);
			hists[h.chip].fillRotation(h, alignment.chips[h.chip].getShiftedCOM() );
			quad.fillHit(h);
			quad.fillRotation(h, alignment.quad.getShiftedCOM());
			auto localPosition=alignment.chips[h.chip].rotateAndShiftBack(h.position);
			auto localResidual=alignment.chips[h.chip].rotateBack(h.residual, {0,0,0} );
			hists[h.chip].local.fill( localPosition, localResidual, h.ToT );
			quad.local.fill( localPosition, localResidual, h.ToT );

			averageHitPosition=averageHitPosition+h.position;
		}
		if(nHitsPassedTotal) averageHitPosition=1./nHitsPassedTotal*averageHitPosition;

		for(auto& h : hists) {
			h.fillEvent();
		}
		quad.fillEvent();

		fitResults.Fill();

		bool draw=true;
		if(draw) {
			std::cout<<"nHitsPassed="<<nHitsPassedTotal<<"";
//			SimpleDetectorConfiguration setupForDrawing { 10,40 /*x*/, 0,42 /*y beam*/, -10,40/*z drift*/};
			auto setupForDrawing=simpleDetectorFromChipCorners(alignment.getAllChipCorners());
			setupForDrawing.minz=-5E3*alignment.driftSpeed.value, setupForDrawing.maxz=5E3*alignment.driftSpeed.value;

//			HoughTransformer::drawCluster(posHits,setupForDrawing);
//			drawQuadOutline(alignment, setupForDrawing.zmin() );

			drawCluster2D(posHits,setupForDrawing);
			alignment.drawChipEdges();

			gPad->Update();
			if(processDrawSignals()) break;
		}
	}
	file.Write();
}
