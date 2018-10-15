/*
 * LaserDataFitter.cpp
 *
 *  Created on: Jul 16, 2018
 *      Author: cligtenb
 */

#include "QuadTrackFitter.h"

#include <algorithm>

#include "TH1.h"
#include "TH2.h"
#include "TTree.h"

#include "PositionHit.h"
#include "../Alignment/TimeWalkCorrector.h"
#include "../TelescopeTrackFitter/HoughTransformer.h"

#pragma link C++ class std::vector<int>+;

namespace {

	const int nChips=4;
//	struct shift {double x,y;};
//	std::array<shift,nChips> chipShifts={{ {24.1,15.7}, {38.0, 26.3}, {23.7,26.3}, {9.9,15.7} }};


	std::vector<PositionHit> convertHitsQuad( TTreeReaderValue<std::vector<Hit> >* chips, double driftSpeed) {
		std::vector<PositionHit> hits;
		for(int i=0; i<nChips; i++) {
			auto chipHits=convertHitsTPC(*chips[i], i, driftSpeed);
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

	struct ChipHistogrammer {
		ChipHistogrammer(const char* name) : dirName(name) {
			using namespace std;

			//change directory
			auto startDir=gDirectory;
			gDirectory->mkdir( dirName.c_str() )->cd();

			nHits=		unique_ptr<TH1I>(new TH1I{"nHits", "Number of hits per event;Number of hits;Entries", 100,0,100});
			constexpr double ToABinWidth=1.5625E-3;
			driftTime=unique_ptr<TH1D>(new TH1D{"driftTime","Drift time;Drift time [#mus];Hits", int(0.8/ToABinWidth),-0.39999,0.4}); //zero would be rounded up (wrong), so move bin slightly to 0.099999
			ToT=			unique_ptr<TH1D>(new TH1D{"ToT", "Time over threshold;ToT [#mus];Hits",80,0,2});
			zHit=			unique_ptr<TH1D>(new TH1D{"zHit", "z-position of hits;z [mm]; Hits", 200,-10,40 });

			xResidual=unique_ptr<TH1D>(new TH1D{"xResidual", "x-residual;x-residual [mm]; Hits", 80,-2,2});
			yResidual=unique_ptr<TH1D>(new TH1D{"yResidual", "y-residual;y-residual [mm]; Hits", 80,-2,2});
			zResidual=unique_ptr<TH1D>(new TH1D{"zResidual", "z-residual;z-residual [mm]; Hits", 80,-4,4});

			pixelHitMap=unique_ptr<TH2D>(new TH2D{"pixelHitMap", "Hitmap by pixel;Columns;Rows", 256,0, 256, 256,0,256});
			positionHitMap=unique_ptr<TH2D>(new TH2D{"positionHitMap", "Hitmap with positions;x-position [mm];y-position [mm]",int(30/.055),11,41,int(40/.055),1,41});

			xRotation=unique_ptr<TH1D>(new TH1D{"xRotation", "x-rotation;x-rotation [rad.]; Hits", 100,-0.5,0.5});
			yRotation=unique_ptr<TH1D>(new TH1D{"yRotation", "y-rotation;y-rotation [rad.]; Hits", 100,-0.5,0.5});
			zRotation=unique_ptr<TH1D>(new TH1D{"zRotation", "z-rotation;z-rotation [rad.]; Hits", 100,-1,1});

		  zResidualByToT=std::unique_ptr<TH2D>(new TH2D{"zResidualByToT", "z-residual by ToT;ToT [#mus]; z-residual [mm]", 100,0,2.5, 200,-5,5});
		  zResidualByToTCorrected=std::unique_ptr<TH2D>(new TH2D{"zResidualByToTCorrected", "z-residual by ToT;ToT [#mus]; z-residual [mm]", 100,0,2.5, 200,-5,5});

		  xResidualByz=std::unique_ptr<TH2D>(new TH2D{"xResidualByz", "x-residuals as a function of drift distance;Drift distance [mm];x-residual [mm]", 26,-1,25,25,-2,2});
		  yResidualByz=std::unique_ptr<TH2D>(new TH2D{"yResidualByz", "y-residuals as a function of drift distance;Drift distance [mm];y-residual [mm]", 26,-1,25,25,-2,2});
		  zResidualByz=std::unique_ptr<TH2D>(new TH2D{"zResidualByz", "z-residuals as a function of drift distance;Drift distance [mm];z-residual [mm]", 26,-1,25,50,-2,2});

			startDir->cd();
		}

		std::unique_ptr<TH1I> nHits;
		std::unique_ptr<TH1D> driftTime, ToT, zHit;
		std::unique_ptr<TH1D> xResidual, yResidual, zResidual;
		std::unique_ptr<TH2D> pixelHitMap, positionHitMap;
		std::unique_ptr<TH1D> xRotation, yRotation, zRotation;
		std::unique_ptr<TH2D> zResidualByToT, zResidualByToTCorrected;
		std::unique_ptr<TH2D> xResidualByz, yResidualByz,zResidualByz;

		void fillHit(const PositionHit& h);
		void fillRotation(const PositionHit& h, const TVector3& COM);
		void fillEvent();
		void fillTimewalk(const PositionHit& h, const TimeWalkCorrector& twc) {
			zResidualByToT->Fill(h.ToT*25E-3, h.residual.z+twc.getCorrection(h.ToT));
		}
	private:
		std::string dirName;
		int hitsCounter=0;
	};




	void ChipHistogrammer::fillHit(const PositionHit& h) {
		driftTime->Fill(h.driftTime/4096.*25E-3);
		ToT->Fill(h.ToT*25E-3);
		zHit->Fill(h.position.z);
		xResidual->Fill(h.residual.x);
		yResidual->Fill(h.residual.y);
		zResidual->Fill(h.residual.z);
		pixelHitMap->Fill(h.column, h.row);
		positionHitMap->Fill(h.position.x, h.position.y);
		zResidualByToTCorrected->Fill(h.ToT*25E-3, h.residual.z);
		xResidualByz->Fill(h.position.z,h.residual.x);
		yResidualByz->Fill(h.position.z,h.residual.y);
		zResidualByz->Fill(h.position.z,h.residual.z);
		hitsCounter++;
	}

	void ChipHistogrammer::fillRotation(const PositionHit& h, const TVector3& COM) {
		//todo: do proper 3d rotation alignment
		auto d=h.position-COM;
		//x
		double perp2x=d.Perp2({1,0,0});
		double phix=(d.z()*h.residual.y-d.y()*h.residual.z)/perp2x;
		xRotation->Fill(phix, perp2x);

		//y
		double perp2y=d.Perp2({0,1,0});
		double phiy=(d.x()*h.residual.z-d.z()*h.residual.x)/perp2y;
		yRotation->Fill(phiy, perp2y);

		//z
		double phi=(d.y()*h.residual.x-d.x()*h.residual.y)/d.Perp2();
		double weight=d.Perp2();
		zRotation->Fill(phi,weight);
	}
	void ChipHistogrammer::fillEvent() {
		if(!hitsCounter) return;
		nHits->Fill(hitsCounter);
		hitsCounter=0;
	}

}



QuadTrackFitter::QuadTrackFitter(std::string fileName) :
		tree("data", fileName) {}

QuadTrackFitter::~QuadTrackFitter() {
}

std::vector<PositionHit> QuadTrackFitter::getSpaceHits(const Alignment& alignment) {
	//retrieve hits
	const auto driftSpeed = alignment.driftSpeed.value / 4096 * 25; //in units of mm/ticks
	posHits = convertHitsQuad(tree.chip, driftSpeed);
	for (auto& h : posHits) {
		//apply alignment
		h = alignment.transform(h);
		h = alignment.timeWalk.correct(h);
		//			if(h.position.z<alignment.hitErrors.z0) h.flag=PositionHit::Flag::smallz;
		h.error = alignment.hitErrors.hitError(h.position.z);
		/*
		 h.calculateResidual(laser);
		 h=flagResidual(h, {2,2,2});
		 h=flagResidualPull(h, {2,2,3});
		 */
	}
	return posHits;
}

void QuadTrackFitter::Loop(std::string outputFile,const Alignment& alignment) {

	TFile file(outputFile.c_str(), "RECREATE");

	Vec3 averageHitPosition;
	std::vector<int> nHits(4);
	int nHitsPassedTotal{0};

	TTree fitResults("fitResults", "fitResults");
	fitResults.Branch("hits", "std::vector<PositionHit>", &posHits);
	fitResults.Branch("hitAverage", "Vec3", &averageHitPosition);
	fitResults.Branch("nHitsPassedChip", "std::vector<int>", &nHits);
	fitResults.Branch("nHitsPassed", &nHitsPassedTotal);

  auto zResidualByToT=std::unique_ptr<TH2D>(new TH2D{"zResidualByToT", "z-residual by ToT;ToT [#mus]; z-residual [mm]", 100,0,2.5, 200,-4,6});
  auto zResidualByToTCorrected=std::unique_ptr<TH2D>(new TH2D{"zResidualByToTCorrected", "z-residual by ToT;ToT [#mus]; z-residual [mm]", 100,0,2.5, 200,-4,6});

	std::array<ChipHistogrammer, nChips> hists{"chip1","chip2","chip3","chip4"};
	ChipHistogrammer quad{"quad"};

	auto nEntries=tree.reader.GetEntries(false);
	tree.reader.SetEntriesRange(1E5,2E5);
	while(tree.reader.Next()) {
		if( !(tree.reader.GetCurrentEntry()%10000) ) std::cout<<"entry "<<tree.reader.GetCurrentEntry()<<" / "<<nEntries<<"\n";
		//retrieve hits
		posHits=getSpaceHits(alignment);
		nHits=std::vector<int>(4);
		for(auto& h : posHits) {
			if(h.flag==PositionHit::Flag::valid) nHits[h.chip]++;
		}
		nHitsPassedTotal=std::accumulate(nHits.begin(), nHits.end(),0);

		averageHitPosition={0,0,0};
		for(auto& h : posHits) {

			if(h.flag==PositionHit::Flag::highResidualxy) continue;
			hists[h.chip].fillTimewalk(h, alignment.timeWalk);
			quad.fillTimewalk(h, alignment.timeWalk);
			zResidualByToTCorrected->Fill(h.ToT*25E-3,h.residual.z);
			zResidualByToT->Fill(h.ToT*25E-3,h.residual.z+alignment.timeWalk.getCorrection(h.ToT));

			if(h.flag!=PositionHit::Flag::valid) continue;
			hists[h.chip].fillHit(h);
			hists[h.chip].fillRotation(h, alignment.chips[h.chip].COM);
			quad.fillHit(h);
			quad.fillRotation(h, alignment.quad.COM);

			averageHitPosition=averageHitPosition+h.position;
		}
		if(nHitsPassedTotal) averageHitPosition=1./nHitsPassedTotal*averageHitPosition;

		for(auto& h : hists) {
			h.fillEvent();
		}
		quad.fillEvent();

		fitResults.Fill();

		SimpleDetectorConfiguration setupForDrawing { 10,40 /*x*/, 0,42 /*y beam*/, 0,40/*z drift*/};
		HoughTransformer::drawCluster(posHits,setupForDrawing);
		if(processDrawSignals()) break;

	}
	file.Write();
}
