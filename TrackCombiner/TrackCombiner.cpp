/*
 * TrackCombiner.cpp
 *
 *  Created on: Oct 15, 2018
 *      Author: cligtenb
 */

#include "TrackCombiner.h"

#include <memory>

#include "deformationCorrection.h"

namespace {

	PositionHit& calculateResidualTimepix(PositionHit& h, const FitResult3D& telescopeFit ) {
		//todo: make residuals perpendicular
		h.residual.x=h.position.x - telescopeFit.xAt(h.position.y);
		h.residual.z=h.position.z - telescopeFit.yAt(h.position.y);
		h.residual.y=0;
		return h;
	}

	std::vector<PositionHit>& calculateResidualTimepix(std::vector<PositionHit>& quadHits,
				const FitResult3D& fit) {
		for (auto& h : quadHits) {
			h = calculateResidualTimepix(h, fit);
		}
		return quadHits;
	}

	bool matchByAverage(std::vector<PositionHit>& quadHits,
			const FitResult3D& fit) {
		//matchin by simple average
		auto avPosition = getAveragePosition(quadHits);
		bool matched=false;
		if (std::fabs(avPosition.x() - fit.xAt(avPosition.y())) < 2
				&& std::fabs(avPosition.z() - fit.yAt(avPosition.y())) < 2) {
			matched = true;
		}
		return matched;
	}

}

TrackCombiner::TrackCombiner(std::string quadFile, std::string telescopeFile, Alignment& alignment) :
	quadFitter(quadFile),
	telescopeFitter(telescopeFile, mimosa),
	alignment(alignment)
{
	//init telescope
	telescopeFitter.setAlignment(alignment);
	telescopeFitter.makeMask(1E4);
	telescopeFitter.maxResidual=0.05;

	//init quad fitter
}

TrackCombiner::~TrackCombiner() {
	outputFile->Write();
}

void TrackCombiner::openFile(std::string outputFileName) {
	outputFile=std::unique_ptr<TFile>(new TFile(outputFileName.c_str(), "RECREATE"));

	const bool makeOutputTree=false;
	if(makeOutputTree) {
		outputTree=std::unique_ptr<TTree>(new TTree("fitResults", "fitResults"));

		outputTree->Branch("telescopeFits", "std::vector<FitResult3D>", &currentEntry.telescopeFits);
		outputTree->Branch("meanQuadPosition", "Vec3", &currentEntry.meanQuadPosition);
		outputTree->Branch("nHitsPerChip", "std::vector<int>", &currentEntry.nHitsPerChip);
		outputTree->Branch("quadHits", "std::vector<PositionHit>", &currentEntry.quadHits);
		outputTree->Branch("matched", &currentEntry.matched);

		outputTree->Branch("telescopeTime", &currentEntry.telescopeTime);
		outputTree->Branch("triggerToA", &currentEntry.triggerToA);
	}


	for(int i=0; i<nChips; i++) {
		hists[i]=std::unique_ptr<ChipHistogrammer>( new ChipHistogrammer("chip"+std::to_string(i), alignment) );
	}
	quadHist=std::unique_ptr<ChipHistogrammer>( new ChipHistogrammer("quad", alignment) );

	auto top=gDirectory;
	top->mkdir("eventCuts")->cd();
	selectedHitAverageToTrackx=std::unique_ptr<TH1D>( new TH1D("selectedHitAverageToTrackx", "distance from average of selected hits to track;Distance [mm];Entries", 100,-1,1) );
	fractionInTrack=std::unique_ptr<TH1D>( new TH1D("fractionInTrack", "fraction of hits in track divided by hits in larger window;Fraction;Entries", 100,0,1) );
	smallestShift=std::unique_ptr<TH1D>( new TH1D("smallestShift", "Smallest shift of all hits in track;Shift [409.6 #mus];",100,0,100) );
	averageShift=std::unique_ptr<TH1D>( new TH1D("averageShift", "Average shift of all hits in track; Shift [409.6 #mus];", 100,0,100 ) );
	top->cd();

}


void TrackCombiner::Process() {
	//match entries
	nTelescopeTriggers=0;
	telescopeFitter.getEntry(previousTriggerNumberBegin);

	quadFitter.reader.tree->GetEntry(0);


	for(int telescopeEntryNumber=0,tpcEntryNumber=0; //5000000, 2308829
			telescopeEntryNumber<2E6 //telescopeFitter.nEvents//1000000
			;) {
//		triggerStatusHistogram.reset();

		// Get Entry and match trigger Numbers
		auto matchStatus=getAndMatchEntries(telescopeEntryNumber,tpcEntryNumber);
		printTriggers(telescopeEntryNumber,tpcEntryNumber);
		if( matchStatus == MatchResult::end ) {
			std::cout<<"\nmatchResult is end\n\n";
			break;}
		else if( matchStatus == MatchResult::noMatch) {
//			std::cout<<" no match\n\n";
			continue; }
//		else { std::cout<<"match!\n\n"; }

//		continue; //debug!!

		auto telescopeHits=telescopeFitter.getSpaceHits();
		int nTelescopeHits=0;
		for(const auto& plane : telescopeHits) nTelescopeHits+=plane.size();
//		std::cout<<"there are "<<nTelescopeHits<<" telescope hits\n";
		if(!nTelescopeHits) continue;

 		auto quadHits=quadFitter.getSpaceHits(alignment);
 		if(!quadHits.size()) continue;
 		auto nHitsPerChip= countHitsPerChip(quadHits);

 		const int minHitsPerChip=10;
 		if( (nHitsPerChip[0]+nHitsPerChip[1] <2*minHitsPerChip) and (nHitsPerChip[2]+ nHitsPerChip[3]<2*minHitsPerChip ) ) continue;

// 		for(auto& h : quadHits) flagShiftedTrigger(h,5); //max shifted

 		//naive fiducial, better below |
// 		for(auto& h : quadHits) {
// 			if(h.row < 32 || h.row >224 || h.column<32 || h.column>224) {
// 				h.flag=PositionHit::Flag::outsideFiducial;
// 			}
// 		}

 		std::vector<FitResult3D> telescopeFits=telescopeFitter.getFits(telescopeHits);
 		std::vector<FitResult3D> timepixFits;

 		bool matched=false;
 		FitResult3D* fittedTrack=nullptr; //todo: save element of array, because if reallocation happens, this pointer is invalidated!
 		for(auto&f : telescopeFits) {

 			//fiducial region
 			const bool useFiducialExpected=false;
 			if(useFiducialExpected) {
				Vec3 expectedAtQuad{ f.xAt(193), 193, f.yAt(193) };
				auto localExpectedAtQuad=alignment.quad.rotateAndShiftBack(expectedAtQuad);
				if( (localExpectedAtQuad.x>11.5 and localExpectedAtQuad.x<17.5) or localExpectedAtQuad.x<2.5 or localExpectedAtQuad.x>25.5 )continue;
 			}

 			//match by calculating number of hits in range
 			const bool matchByResidual=true;
 			if(matchByResidual) {


 	 			auto quadHitsWithResidual=quadHits; //copy hits!
 	 			//cut with fixed window
				for(auto& h : quadHitsWithResidual) {
					h=calculateResidualTimepix(h,f);
					h=flagResidual(h,{1.5,1.5,2});
				}

				int nHitsAfterSimpleCut=countTotalValidHits(quadHitsWithResidual);
				if(nHitsAfterSimpleCut<20) {
					if(nHitsAfterSimpleCut && drawEvent) std::cout<<nHitsAfterSimpleCut<<" is less than 20 hits along track!\n";
					continue;
				}

				//convert hits to telescope coordinates x,y,z -> x,z,y
				std::vector<PositionHit> quadHitsInTelescopeCoords;
				for(const auto& h : quadHitsWithResidual) {
					PositionHit th{h};
					std::swap(th.position.y, th.position.z);
					std::swap(th.error.y, th.error.z);
					quadHitsInTelescopeCoords.push_back(th);
				}
				//fit hits
 	 			FitResult3D timepixFit=regressionFit3d(quadHitsInTelescopeCoords,0);
 	 			timepixFits.push_back(timepixFit);
 	 			for(auto& h : quadHitsWithResidual) {
 	 				Vec3 timepixExpectedPos{h.position.x, h.position.y, timepixFit.yAt(h.position.z)};
 	 				h.error=alignment.hitErrors.hitError( alignment.quad.rotateAndShiftBack(timepixExpectedPos).z );
					h=flagResidualPull(h, {2,5,3} );
 	 			}


				int nTotalValidHits=countTotalValidHits(quadHitsWithResidual);
				if(nTotalValidHits<20) {
					if(nTotalValidHits && drawEvent) std::cout<<nTotalValidHits<<" is less than 20 hits along track!\n";
					continue;
				}

				//fraction
				int nSideBandHits=std::count_if(quadHitsWithResidual.begin(), quadHitsWithResidual.end(), [](const PositionHit& h) {
					return h.flag==PositionHit::Flag::valid
							or ( (h.flag==PositionHit::Flag::highResidualxy or h.flag==PositionHit::Flag::highResidualz)
//									and fabs(h.residual.x/h.error.x)<4 and fabs(h.residual.z/h.error.z)<4);
									and fabs(h.residual.x)<5 and fabs(h.residual.z)<2);
				});
				double fraction=double(nHitsAfterSimpleCut)/nSideBandHits;
				fractionInTrack->Fill(fraction);

				//average pos
				auto avPosition = getAveragePosition(quadHitsWithResidual,true);
				double averageDist=avPosition.x() - f.xAt(avPosition.y()) ;
				selectedHitAverageToTrackx->Fill(averageDist);

				//first occurence (shift) and average shift
				if(drawEvent) std::cout<<"dist "<<averageDist<<" frac "<<fraction<<"\n";
				if ( std::fabs(averageDist) > 1 || fraction<0.8 ) {
					if(drawEvent) std::cout<<"More than 20 hits along track, but average position did not match or fraction to low!\n";
					continue;
				}


				int sumShift=0, nValid=0;
				short minShift=quadHitsWithResidual.front().nShiftedTrigger;
				{	for(const auto& h : quadHitsWithResidual) {
						if(!h.isValid()) continue;
						sumShift+=h.nShiftedTrigger;
						minShift=std::min(minShift, h.nShiftedTrigger);
						nValid++;
					}
					smallestShift->Fill(minShift);
					averageShift->Fill(sumShift/nValid);
				}

				if ( minShift>5 || sumShift/nValid >150 ) {
					if(drawEvent) std::cout<<"More than 20 hits along track, but first hit shifted by too much or average shift is too high!\n";
					continue;
				}


				matched = true;
				fittedTrack=&f;
				quadHits=quadHitsWithResidual;
				break;

 			} else {
 			//matchin by simple average
				if( matchByAverage(quadHits, f) ) {
					matched=true;
					quadHits=calculateResidualTimepix(quadHits,f);
					fittedTrack=&f;
					break;
				}
 			}
 		}
// 		if(!matched) continue;


 		if(matched) {
 			for(auto& h : quadHits) {
 				if(h.flag!=PositionHit::Flag::valid) continue;
 				hists[h.chip]->fillHit(h);
 				hists[h.chip]->fillRotation(h, alignment.chips[h.chip].getCOMGlobal());
 				quadHist->fillHit(h);
 				quadHist->fillRotation(h, alignment.quad.getCOMGlobal());

 				auto localPosition=alignment.quad.rotateAndShiftBack(h.position);
 				auto expectedPosition = Vec3{ fittedTrack->xAt(h.position.y), h.position.y, fittedTrack->yAt(h.position.y)}; //telescope->timepix frame!
 				auto localExpectedPosition=alignment.quad.rotateAndShiftBack(expectedPosition);
 				auto localResidual=alignment.quad.rotateBack(h.residual, {0,0,0});

 				hists[h.chip]->fillTimewalk(localResidual.z, h.ToT, alignment.timeWalk);
 				quadHist->fillTimewalk(localResidual.z, h.ToT, alignment.timeWalk);

 				hists[h.chip]->local.fill( localPosition, localResidual, h.ToT);
 				quadHist->local.fill( localPosition, localResidual, h.ToT );
 				quadHist->locExp.fill(localExpectedPosition, localResidual, h.ToT);
 				hists[h.chip]->locExp.fill( localExpectedPosition, localResidual, h.ToT );
 				if(fabs(h.residual.z/h.error.z)<1.5) {
 					hists[h.chip]->tighterCuts.fill( localExpectedPosition, localResidual, h.ToT );
 					quadHist->tighterCuts.fill( localExpectedPosition, localResidual, h.ToT );
 				}

 				auto correctedResidual=localResidual;
 				correctedResidual.x-=deformationCorrection( h.chip, localExpectedPosition.x );//using _expected_ position.
 				hists[h.chip]->corrected.fill( localExpectedPosition, correctedResidual, h.ToT );
 				quadHist->corrected.fill(localExpectedPosition, correctedResidual, h.ToT);
 			}
 		}


		for(auto& h : hists) {
			h->fillEvent();
		}
		quadHist->fillEvent();

		if(drawEvent) {
			//		SimpleDetectorConfiguration setupForDrawing { 0,30 /*x*/, 0,42 /*y beam*/, -20,20/*z drift*/};
			auto setupForDrawing=simpleDetectorFromChipCorners(alignment.getAllChipCorners());
			setupForDrawing.minz=-10, setupForDrawing.maxz=30;
			//		SimpleDetectorConfiguration setupForDrawing { 10,40 /*x*/, 0,400 /*y beam*/, 0,40/*z drift*/};

//			std::cout<<"draw telescope hits\n";
//			telescopeFitter.drawEvent(telescopeHits,telescopeFits);
	//		HoughTransformer::drawClusters(telescopeHits, setupForDrawing);

			std::cout<<"draw "<< (matched?"matched":"") <<" quad hits "<<quadHits.size()<<"\n";

			//3D
			const bool draw3D=true;
			if(draw3D) {
				HoughTransformer::drawCluster(quadHits,setupForDrawing);
				for (auto& f : telescopeFits)
					f.draw( setupForDrawing.ymin(), setupForDrawing.ymax() );
				for (auto& f : timepixFits)
					f.draw( setupForDrawing.ymin(), setupForDrawing.ymax(), kTeal );
				drawQuadOutline(alignment, setupForDrawing.zmax() );
				gPad->Update();
			}
			const bool draw2D=true;
			if(draw2D){
				//2D
				drawCluster2D(quadHits,setupForDrawing);
				alignment.drawChipEdges();
				for (auto& f : telescopeFits)
					f.XZ.draw( setupForDrawing.ymin(), setupForDrawing.ymax() );
				gPad->Update();
			}


			gPad->Update();
			if(processDrawSignals()) break;
		}

		if(outputTree) {
			currentEntry.telescopeFits=telescopeFits;
			currentEntry.meanQuadPosition=getAveragePosition(quadHits, {PositionHit::Flag::shiftedTrigger});
			currentEntry.nHitsPerChip=nHitsPerChip;
			currentEntry.quadHits=quadHits;
			currentEntry.matched=matched;

			currentEntry.telescopeTime=telescopeFitter.timestamp;
			currentEntry.triggerToA=quadFitter.reader.triggerToA;

			outputTree->Fill();
		}

	}
	std::cout<<"\ndone!\n";

	if(not drawEvent and outputTree) {
		static TCanvas* xCorrelation= new TCanvas();
		xCorrelation->cd();
		outputTree->Draw("XZ.intercept+XZ.slope*172:x", "", "", 5E4);
		static TCanvas* yCorrelation= new TCanvas();
		yCorrelation->cd();
		outputTree->Draw("YZ.intercept+YZ.slope*172:z", "fabs(XZ.intercept+XZ.slope*172-x)<3", "", 5E4);
	}
}

void TrackCombiner::printTriggers(int telescopeEntry, int tpcEntry) {
	int printEveryN=1000;
	static int i=0;
	using namespace std;
//	if(telescopeEntry>287035) printEveryN=1;
	if( !(i++%printEveryN) ) {
		cout<<"entry: "<<telescopeEntry<<"/"<<telescopeFitter.nEvents<<" ";
		cout<<"triggers: "<<telescopeFitter.triggerNumberBegin<<"-"<<telescopeFitter.triggerNumberEnd;
		cout<<" timepix triggerNumber: "<<quadFitter.triggerNumber()<<"="<<(quadFitter.triggerNumber()+triggerOffset) % 32768<<" in entry "<<tpcEntry<<"\r"<<flush;
	}
}

TrackCombiner::MatchResult TrackCombiner::getAndMatchEntries(
		int& telescopeEntry,
		int& tpcStartEntry) {

	static bool printMatching = false;
//	if(telescopeEntry>287035) printMatching=true;
//	if(telescopeEntry>805) printMatching=false;

	if(printMatching) std::cout<<"starting getAndMatchEntries( "<<telescopeEntry<<", "<<tpcStartEntry<<")\n";

	//use modulus offset instead of doing actual modulus, to continue counting after 32768
	long long modulusOffset=32768*int((quadFitter.triggerNumber()+triggerOffset)/32768);

	//special case: triggerNumberBegin==triggerNumberEnd==0
	while(telescopeFitter.triggerNumberBegin==0 and telescopeFitter.triggerNumberEnd==0) {
		std::cout<<"special case telescope trigger number begin==end==0\n";
		if( !telescopeFitter.getEntry(++telescopeEntry) ) return MatchResult::end;
	}

	//if telescope triggerNumberBegin decreased, and tpc did not already pass this boundary
	//then we must get entries until the tpc triggernumber also passes the 32768 boundary ( is larger or equal to new trigger number)
	if(telescopeFitter.triggerNumberBegin<previousTriggerNumberBegin-10
			and not (quadFitter.triggerNumber()+triggerOffset-modulusOffset<500) ) {
			if(printMatching)		std::cout<<"telescope number decreased and tpc did not yet, get entries until tpc also passes boundary\n";
		do {
			if(printMatching)			printTriggers(telescopeEntry, tpcStartEntry);
			if( not quadFitter.getEntry(tpcStartEntry++) ) return MatchResult::end;
		} while ( quadFitter.triggerNumber()+triggerOffset - modulusOffset - 32768 < telescopeFitter.triggerNumberBegin );
		modulusOffset=32768*int((quadFitter.triggerNumber()+triggerOffset)/32768);
		//stop timepix from going back over this boundary
		timepixEntryFirstMatch=tpcStartEntry;
		hadFirstMatch=true;
	}

	//get next entry until tpc trigger number is larger than or equal to begin
	//or until the tpc trigger number decreases;
	do {
		if(printMatching)		std::cout<<"getting tpc, until tpc trigger entry is larger or equal to begin\n";
		if( not quadFitter.getEntry(tpcStartEntry++) ) {
			if(printMatching) std::cout<<"could not get entry "<<tpcStartEntry-1<<"\n";
			return MatchResult::end;
		}
	} while( (quadFitter.triggerNumber()+triggerOffset-modulusOffset) < telescopeFitter.triggerNumberBegin
			and telescopeFitter.triggerNumberBegin - (quadFitter.triggerNumber()+triggerOffset-modulusOffset)<500 );

	//if also larger than end: reached the end of this telescope frame, continue with next telescope frame;
	if( (quadFitter.triggerNumber()+triggerOffset-modulusOffset) > telescopeFitter.triggerNumberEnd) {
		if(printMatching)		std::cout<<"continuing to next telescope frame\n";
//		triggerStatusHistogram.Fill("Trigger numbers do not match", 1);

//		frameStatusHistogram.reset();
		nTelescopeTriggers+=std::max(0,telescopeFitter.triggerNumberEnd-previousTriggerNumberEnd);

		previousTriggerNumberBegin=telescopeFitter.triggerNumberBegin;
		previousTriggerNumberEnd=telescopeFitter.triggerNumberEnd;
		if( !telescopeFitter.getEntry(++telescopeEntry) ) return MatchResult::end;

		if(printMatching)		std::cout<<"increased telescopeEntry, first time pix match was: "<<timepixEntryFirstMatch<<std::endl;

		tpcStartEntry=timepixEntryFirstMatch;
		hadFirstMatch=false;

		return MatchResult::noMatch;
	}

	if(not hadFirstMatch) {
		timepixEntryFirstMatch=tpcStartEntry-1;
		hadFirstMatch=true;
	}

	return MatchResult::match;

}
