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

	std::vector<PositionHit> getInTelescopeCoords( std::vector<PositionHit> hits ) {
		std::vector<PositionHit> quadHitsInTelescopeCoords;
		for(const auto& h : hits) {
			PositionHit th{h};
			std::swap(th.position.y, th.position.z);
			std::swap(th.error.y, th.error.z);
			quadHitsInTelescopeCoords.push_back(th);
		}
		return quadHitsInTelescopeCoords;
	}

}

TrackCombiner::TrackCombiner(std::string quadFile, std::string telescopeFile, Alignment& alignment) :
	quadFitter(quadFile),
	telescopeFitter(telescopeFile, mimosa),
	alignment(alignment)
{
	//init telescope
	telescopeFitter.setAlignment(alignment);
	telescopeFitter.makeMask(1E5);
	telescopeFitter.maxResidual=0.015;

	//init quad fitter
}

TrackCombiner::~TrackCombiner() {
	outputFile->Write();
}

void TrackCombiner::openFile(std::string outputFileName) {
	outputFile=std::unique_ptr<TFile>(new TFile(outputFileName.c_str(), "RECREATE"));

	const bool makeOutputTree=true;
	if(makeOutputTree) {
		outputTree=std::unique_ptr<TTree>(new TTree("fitResults", "fitResults"));

		outputTree->Branch("telescopeFits", "std::vector<FitResult3D>", &currentEntry.telescopeFits);
		outputTree->Branch("timepixFits", "std::vector<FitResult3D>", &currentEntry.timepixFits);
		outputTree->Branch("meanQuadPosition", "Vec3", &currentEntry.meanQuadPosition);
		outputTree->Branch("meanPositionPerChip", "std::vector<Vec3>", &currentEntry.meanPositionPerChip);
		outputTree->Branch("meanQuadDiff", "Vec3", &currentEntry.meanQuadDiff);
		outputTree->Branch("meanDiffPerChip", "std::vector<Vec3>", &currentEntry.meanDiffPerChip);
		outputTree->Branch("meanDiffPerChipPerFitFirst", "std::vector<std::vector<Vec3>>", &currentEntry.meanDiffPerChipPerFitFirst);
		outputTree->Branch("meanDiffPerChipPerFitLast", "std::vector<std::vector<Vec3>>", &currentEntry.meanDiffPerChipPerFitLast);
		outputTree->Branch("centerDiffPerChipFit", "Vec3", &currentEntry.centerDiffPerChipFit);
		outputTree->Branch("centerErrorPerChipFit", "Vec3", &currentEntry.centerErrorPerChipFit);
		outputTree->Branch("meanErrorPerChipPerFitFirst", "std::vector<std::vector<Vec3>>", &currentEntry.meanErrorPerChipPerFitFirst);
		outputTree->Branch("meanErrorPerChipPerFitLast", "std::vector<std::vector<Vec3>>", &currentEntry.meanErrorPerChipPerFitLast);
		outputTree->Branch("meanQuadError", "Vec3", &currentEntry.meanQuadError);
		outputTree->Branch("meanErrorPerChip", "std::vector<Vec3>", &currentEntry.meanErrorPerChip);
		outputTree->Branch("nHitsPerChip", "std::vector<int>", &currentEntry.nHitsPerChip);
		outputTree->Branch("nHitsPerChipValid", "std::vector<int>", &currentEntry.nHitsPerChipValid);
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
	selectedHitAverageToTrackz=std::unique_ptr<TH1D>( new TH1D("selectedHitAverageToTrackz", "distance from average of selected hits to track;Distance [mm];Entries", 100,-1,1) );
	fractionInTrack=std::unique_ptr<TH1D>( new TH1D("fractionInTrack", "fraction of hits in track divided by hits in larger window;Fraction;Entries", 100,0,1) );
	smallestShift=std::unique_ptr<TH1D>( new TH1D("smallestShift", "Smallest shift of all hits in track;Shift [409.6 #mus];",100,0,100) );
	averageShift=std::unique_ptr<TH1D>( new TH1D("averageShift", "Average shift of all hits in track; Shift [409.6 #mus];", 100,0,100 ) );

	if(keepStatus) {
		frameStatusHistogram=unique_ptr<StatusKeeper>{new StatusKeeper("frame")};
		triggerStatusHistogram=unique_ptr<StatusKeeper>{new StatusKeeper("trigger")};
		timepixStatusHistogram=unique_ptr<StatusKeeper>{new StatusKeeper("timepixTrigger")};
	}

	top->cd();
	top->mkdir("eventhist")->cd();
	averageToTrackx=std::unique_ptr<TH1D>( new TH1D("averageToTrackx", "distance from average of selected hits to track;Distance [mm];Entries", 100,-1,1) );
	top->cd();
}


void TrackCombiner::Process() {
	//match entries
	nTelescopeTriggers=0;
	telescopeFitter.getEntry(previousTriggerNumberBegin);

	quadFitter.reader.tree->GetEntry(0);

	std::cout<<"number of events: telescope="<<telescopeFitter.nEvents<<", timepix="<<quadFitter.numberOfEntries()<<"\n";

	bool pdfIsOpen=false;

	for(int telescopeEntryNumber=0,tpcEntryNumber=0; //5000000, 2308829
			telescopeEntryNumber<1E6 //telescopeFitter.nEvents//1000000
			;) {
		if(keepStatus) triggerStatusHistogram->reset();

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


		auto telescopeHits=telescopeFitter.getSpaceHits();
		int nTelescopeHits=0;
		for(const auto& plane : telescopeHits) nTelescopeHits+=plane.size();
//		std::cout<<"there are "<<nTelescopeHits<<" telescope hits\n";
		if(!nTelescopeHits) {
 			replaceStatus(1,"no telescope hits", tpcEntryNumber);
			continue;
		}

 		auto quadHits=quadFitter.getSpaceHits(alignment);
 		if(!quadHits.size()) {
 			replaceStatus(2,"no quad hits", tpcEntryNumber);
 			continue;
 		}
 		auto nHitsPerChip= countHitsPerChip(quadHits);

 		const int minHitsPerChip=10;
 		if( (nHitsPerChip[0]+nHitsPerChip[1] <2*minHitsPerChip) and (nHitsPerChip[2]+ nHitsPerChip[3]<2*minHitsPerChip ) ) {
 			replaceStatus(10,"Less than 20 hits on a chip pair", tpcEntryNumber);
 			continue;
 		}


// 		for(auto& h : quadHits) flagShiftedTrigger(h,5); //max shifted

 		//naive fiducial, better below |
// 		for(auto& h : quadHits) {
// 			if(h.row < 32 || h.row >224 || h.column<32 || h.column>224) {
// 				h.flag=PositionHit::Flag::outsideFiducial;
// 			}
// 		}

 		std::vector<FitResult3D> telescopeFits=telescopeFitter.getFits(telescopeHits);
 		auto partialTelescopeFits=telescopeFitter.partialFits;
 		std::vector<FitResult3D> timepixFits;

		if(telescopeFits.empty()) { replaceStatus(11, "All telescope clusters failed fit", tpcEntryNumber); continue; }

 		bool matched=false;
 		FitResult3D *fittedTrack{nullptr}, timepixTrack{}, *firstPartialFit{nullptr}, *secondPartialFit{nullptr}; //todo: save element of array, because if reallocation happens, this pointer is invalidated!
 		for(unsigned iTelescopeFit=0; iTelescopeFit<telescopeFits.size();iTelescopeFit++) {
 			auto& f = telescopeFits[iTelescopeFit];

 			//fiducial region
 			const bool useFiducialExpected=true;
 			if(useFiducialExpected) {
				Vec3 expectedAtQuad{ f.xAt(193), 193, f.yAt(193) };
				auto localExpectedAtQuad=alignment.quad.rotateAndShiftBack(expectedAtQuad);
				if( (localExpectedAtQuad.x>11.5 and localExpectedAtQuad.x<17.5) or localExpectedAtQuad.x<2.5 or localExpectedAtQuad.x>25.5 ) {
					replaceStatus(21, "Track does not passs through the fiducial region", tpcEntryNumber);
					continue;
				}
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
					replaceStatus(22, "Less than 20 hits after plain residual cut", tpcEntryNumber);
					continue;
				}

				//convert hits to telescope coordinates x,y,z -> x,z,y
				auto quadHitsInTelescopeCoords=getInTelescopeCoords(quadHitsWithResidual);
				//fit hits
 	 			FitResult3D timepixFit=regressionFit3d(quadHitsInTelescopeCoords);
 	 			timepixFits.push_back(timepixFit);
 	 			//set error and calculate resolution
 	 			for(auto& h : quadHitsWithResidual) {
// 	 				Vec3 timepixExpectedPos{h.position.x, h.position.y, timepixFit.yAt(h.position.z)}; //use timepix fit
 	 				Vec3 timepixExpectedPos{h.position.x, h.position.y, f.yAt(h.position.z)}; //use telescope fit!
 	 				h.error=alignment.hitErrors.hitError( alignment.quad.rotateAndShiftBack(timepixExpectedPos).z );
					h=flagResidualPull(h, {3,5,3} );
 	 			}


				int nTotalValidHits=countTotalValidHits(quadHitsWithResidual);
				if(nTotalValidHits<20) {
					if(nTotalValidHits && drawEvent) std::cout<<nTotalValidHits<<" is less than 20 hits along track!\n";
					replaceStatus(23, "Less than 20 hits after pull residual cut", tpcEntryNumber);
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
				auto avPosition =getWeightedAveragePosition(quadHitsWithResidual, [](const PositionHit& h){return h.isValid();}); //reject all flagged
				double averageDistx=avPosition.x() - f.xAt(avPosition.y()) ;
				double averageDistz=avPosition.z() - f.yAt(avPosition.y()) ;
				selectedHitAverageToTrackz->Fill(averageDistz);

				//first occurence (shift) and average shift
				if(drawEvent) std::cout<<"dist "<<averageDistx<<", "<<averageDistz<<" frac "<<fraction<<"\n";
				if ( std::fabs(averageDistz) > 0.3 || fraction<0.8 ) {
					if(drawEvent) std::cout<<"More than 20 hits along track, but average z position did not match or fraction to low!\n";
					replaceStatus(24, "Telescope track does not match Gridpix hits (y or fraction)", tpcEntryNumber);
					continue;
				}
				selectedHitAverageToTrackx->Fill(averageDistx);
				if ( std::fabs(averageDistx) > 0.3 ) {
					if(drawEvent) std::cout<<"More than 20 hits along track, but average x position did not match!\n";
					replaceStatus(25, "Telescope track does not match Gridpix hits (x)", tpcEntryNumber);
					continue;
				}


				int sumShift=0, nValid=0;
				short minShift=1000;
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
					replaceStatus(26, "Shifted trigger cut", tpcEntryNumber);
					continue;
				}
				if(drawEvent) std::cout<<"min shift is "<<minShift<<" average is "<<sumShift/nValid<<"\n";

				matched = true;
				replaceStatus(100, "Successful", tpcEntryNumber);

				const bool useTimepixFit=false;
				if(useTimepixFit) {
					fittedTrack=&timepixFit;
					quadHitsWithResidual=calculateResidualTimepix(quadHitsWithResidual,timepixFit); //recalculate residuals
				} else {
					fittedTrack=&f;
				}
				firstPartialFit=&partialTelescopeFits[iTelescopeFit].first;
				secondPartialFit=&partialTelescopeFits[iTelescopeFit].second;
				quadHits=quadHitsWithResidual;
				timepixTrack=timepixFit;
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
 				quadHist->fillHit(h);
 				quadHist->fillRotation(h, alignment.quad.getShiftedCOM());

 				auto localPosition=alignment.quad.rotateAndShiftBack(h.position);
 				auto expectedPosition = Vec3{ fittedTrack->xAt(h.position.y), h.position.y, fittedTrack->yAt(h.position.y)}; //telescope->timepix frame!
 				auto localExpectedPosition=alignment.quad.rotateAndShiftBack(expectedPosition);
 				auto localResidual=alignment.quad.rotateBack(h.residual, {0,0,0});

 				auto chipFrameCOM=alignment.chips[h.chip].COM;
 				auto chipFrameResidual=alignment.chips[h.chip].rotateBack(localResidual, {0,0,0});
 				auto chipFrameExpectedPosition=alignment.chips[h.chip].rotateAndShiftBack(localExpectedPosition);
 				hists[h.chip]->fillRotation(chipFrameExpectedPosition, chipFrameResidual, chipFrameCOM );

 				hists[h.chip]->fillTimewalk(localResidual.z, h.ToT, alignment.timeWalks[h.chip]);
 				quadHist->fillTimewalk(localResidual.z, h.ToT, alignment.timeWalks[h.chip]);

 				hists[h.chip]->local.fill( localPosition, localResidual, h.ToT);
 				quadHist->local.fill( localPosition, localResidual, h.ToT );
 				quadHist->locExp.fill(localExpectedPosition, localResidual, h.ToT);
 				hists[h.chip]->locExp.fill( localExpectedPosition, localResidual, h.ToT );
 				if(fabs(h.residual.z/h.error.z)<1.5) {
 					hists[h.chip]->tighterCuts.fill( localExpectedPosition, localResidual, h.ToT );
 					quadHist->tighterCuts.fill( localExpectedPosition, localResidual, h.ToT );
 				} else {
 					hists[h.chip]->sideBands.fill( localExpectedPosition, localResidual, h.ToT );
 					quadHist->sideBands.fill( localExpectedPosition, localResidual, h.ToT );
 				}

 				auto correctedResidual=localResidual;
// 				correctedResidual.x-=deformationCorrection( h.chip, localExpectedPosition.x );//using _expected_ position.
 				correctedResidual.x-=deformationCorrection2D( h.chip, localExpectedPosition.x, localExpectedPosition.y );//using _expected_ position.
 				hists[h.chip]->corrected.fill( localExpectedPosition, correctedResidual, h.ToT );
 				quadHist->corrected.fill(localExpectedPosition, correctedResidual, h.ToT);

 				const bool correctHits=false;
 				if(correctHits) {
 					h.residual=alignment.quad.rotate(localResidual);
 					auto localCorrectedPosition=localPosition;
 					localCorrectedPosition.x-=deformationCorrection2D( h.chip, localExpectedPosition.x, localExpectedPosition.y );
 					h.position=alignment.quad.rotateAndShift(localCorrectedPosition);
 				}
 			}

			auto averagePos=getAveragePosition(quadHits, true);
			averageToTrackx->Fill( averagePos.x()-fittedTrack->xAt(averagePos.y()) );
 		}


		for(auto& h : hists) {
			h->fillEvent();
		}
		quadHist->fillEvent();


		if(keepStatus) {
			int oldestRelevantTpcEntryNumber=tpcEntryNumber-(telescopeFitter.triggerNumberEnd-telescopeFitter.triggerNumberBegin);
			timepixStatusKeepers.writeBufferUpTo(tpcEntryNumber-1E4, [](StatusKeeper&s){s.reset();});
		}

		if(outputTree and matched) {
			//make timepix fit per chip
			auto hitsPerChip=getHitsPerChip(getInTelescopeCoords(quadHits));
			std::vector<std::unique_ptr<FitResult3D>> fitsPerChip, fitsPerChipWithFirstPlanes, fitsPerChipWithLastPlanes;
			//make telescope points
			double firstPlanesz=mimosa.planePosition[2], lastPlanesz=mimosa.planePosition[3];
			PositionHit firstThreePlanes{ firstPartialFit->xAt(firstPlanesz), firstPartialFit->yAt(firstPlanesz), firstPlanesz,0};
			firstThreePlanes.error=Vec3{0.020,0.020,1};
			PositionHit lastThreePlanes{ secondPartialFit->xAt(lastPlanesz), secondPartialFit->yAt(lastPlanesz), lastPlanesz,0};
			lastThreePlanes.error=Vec3{0.020,0.020,1};
			std::vector<PositionHit> firstAndLastIntercept{firstThreePlanes, lastThreePlanes};
			FitResult3D firstToLastPartialFit=regressionFit3d(firstAndLastIntercept);
 			for(auto& hv : hitsPerChip) {
 				if(hv.size()>=2) {
 					fitsPerChip.emplace_back( new FitResult3D(regressionFit3d(hv)) );

 					hv.push_back(firstThreePlanes);
 					fitsPerChipWithFirstPlanes.emplace_back( new FitResult3D(regressionFit3d(hv)) );

 					hv.back()=lastThreePlanes;
 					fitsPerChipWithLastPlanes.emplace_back( new FitResult3D(regressionFit3d(hv)) );

 				} else {
 					fitsPerChip.push_back( nullptr );
 					fitsPerChipWithFirstPlanes.push_back(nullptr);
 					fitsPerChipWithLastPlanes.push_back(nullptr);
 				}
			}

// 			auto fittedTrack=&firstToLastPartialFit; //DEBUG SHADOW actual fittedTrack!!

			currentEntry.telescopeFits=telescopeFits;
			currentEntry.timepixFits=timepixFits;
//			currentEntry.meanQuadPosition=getAveragePosition(quadHits, {PositionHit::Flag::shiftedTrigger});
			currentEntry.meanQuadPosition=getWeightedAveragePosition(quadHits, [](const PositionHit& h){return h.isValid();}); //reject all flagged
			currentEntry.meanPositionPerChip=getWeightedAveragePositionPerChip(quadHits);
			currentEntry.meanQuadDiff.x=matched ? (currentEntry.meanQuadPosition.x - fittedTrack->xAt(currentEntry.meanQuadPosition.y)) : 0;
			currentEntry.meanQuadDiff.z=matched ? (currentEntry.meanQuadPosition.z - fittedTrack->yAt(currentEntry.meanQuadPosition.y)) : 0;
			currentEntry.meanQuadError.x=matched ? hypot(timepixTrack.XZ.errorAt(currentEntry.meanQuadPosition.y), fittedTrack->XZ.errorAt(currentEntry.meanQuadPosition.y)) : 0;
			currentEntry.meanQuadError.z=matched ? hypot(timepixTrack.YZ.errorAt(currentEntry.meanQuadPosition.y), fittedTrack->YZ.errorAt(currentEntry.meanQuadPosition.y)) : 0;
			currentEntry.meanDiffPerChip=std::vector<Vec3>(4);
			currentEntry.meanErrorPerChip=std::vector<Vec3>(4);
			currentEntry.meanDiffPerChipPerFitFirst=std::vector<std::vector<Vec3>>{4,std::vector<Vec3>{4}};
			currentEntry.meanDiffPerChipPerFitLast=std::vector<std::vector<Vec3>>{4,std::vector<Vec3>{4}};
			currentEntry.meanErrorPerChipPerFitFirst=std::vector<std::vector<Vec3>>{4,std::vector<Vec3>{4}};
			currentEntry.meanErrorPerChipPerFitLast=std::vector<std::vector<Vec3>>{4,std::vector<Vec3>{4}};
			for(int i=0; i<4; i++) { //chip with position
				currentEntry.meanDiffPerChip[i].x =	matched ? (currentEntry.meanPositionPerChip[i].x - fittedTrack->xAt(currentEntry.meanPositionPerChip[i].y) ) : 0;
				currentEntry.meanDiffPerChip[i].z =	matched ? (currentEntry.meanPositionPerChip[i].z - fittedTrack->yAt(currentEntry.meanPositionPerChip[i].y) ) : 0;
				currentEntry.meanErrorPerChip[i].x = matched && fitsPerChip[i] ? hypot(fitsPerChip[i]->XZ.errorAt(currentEntry.meanPositionPerChip[i].y),fittedTrack->XZ.errorAt(currentEntry.meanPositionPerChip[i].y) ) : 0;
				currentEntry.meanErrorPerChip[i].z = matched && fitsPerChip[i] ? hypot(fitsPerChip[i]->YZ.errorAt(currentEntry.meanPositionPerChip[i].y),fittedTrack->YZ.errorAt(currentEntry.meanPositionPerChip[i].y)) : 0;
				for(int j=0; j<4; j++) { //chip with fit
					currentEntry.meanDiffPerChipPerFitFirst[i][j].x=fitsPerChipWithFirstPlanes[j] ? (currentEntry.meanPositionPerChip[i].x - fitsPerChipWithFirstPlanes[j]->xAt(currentEntry.meanPositionPerChip[i].y) ) : 0;
					currentEntry.meanDiffPerChipPerFitFirst[i][j].z=fitsPerChipWithFirstPlanes[j] ? (currentEntry.meanPositionPerChip[i].z - fitsPerChipWithFirstPlanes[j]->yAt(currentEntry.meanPositionPerChip[i].y) ) : 0;

					currentEntry.meanDiffPerChipPerFitLast[i][j].x=fitsPerChipWithLastPlanes[j] ? (currentEntry.meanPositionPerChip[i].x - fitsPerChipWithLastPlanes[j]->xAt(currentEntry.meanPositionPerChip[i].y) ) : 0;
					currentEntry.meanDiffPerChipPerFitLast[i][j].z=fitsPerChipWithLastPlanes[j] ? (currentEntry.meanPositionPerChip[i].z - fitsPerChipWithLastPlanes[j]->yAt(currentEntry.meanPositionPerChip[i].y) ) : 0;

					currentEntry.meanErrorPerChipPerFitFirst[i][j].x=fitsPerChipWithFirstPlanes[j] ?
							hypot(fitsPerChipWithFirstPlanes[j]->XZ.errorAt(currentEntry.meanPositionPerChip[i].y),fittedTrack->XZ.errorAt(currentEntry.meanPositionPerChip[i].y)  ) : 0;
					currentEntry.meanErrorPerChipPerFitFirst[i][j].z=fitsPerChipWithFirstPlanes[j] ?
							hypot(fitsPerChipWithFirstPlanes[j]->YZ.errorAt(currentEntry.meanPositionPerChip[i].y),fittedTrack->YZ.errorAt(currentEntry.meanPositionPerChip[i].y) ) : 0;

					currentEntry.meanErrorPerChipPerFitLast[i][j].x=fitsPerChipWithLastPlanes[j] ?
							hypot(fitsPerChipWithLastPlanes[j]->XZ.errorAt(currentEntry.meanPositionPerChip[i].y),fittedTrack->XZ.errorAt(currentEntry.meanPositionPerChip[i].y) ) : 0;
					currentEntry.meanErrorPerChipPerFitLast[i][j].z=fitsPerChipWithLastPlanes[j] ?
							hypot(fitsPerChipWithLastPlanes[j]->YZ.errorAt(currentEntry.meanPositionPerChip[i].y),fittedTrack->YZ.errorAt(currentEntry.meanPositionPerChip[i].y) ) : 0;

				}
			}
			auto quadCenter=alignment.quad.getShiftedCOM().y();
			if(fitsPerChip[0] && fitsPerChip[1]) {
				currentEntry.centerDiffPerChipFit.x = fitsPerChipWithFirstPlanes[1]->xAt(quadCenter) - fitsPerChipWithLastPlanes[0]->xAt(quadCenter);
				currentEntry.centerDiffPerChipFit.z = fitsPerChipWithFirstPlanes[1]->yAt(quadCenter) - fitsPerChipWithLastPlanes[0]->yAt(quadCenter);
				currentEntry.centerErrorPerChipFit.x = sqrt( fitsPerChipWithFirstPlanes[1]->XZ.error2At(quadCenter) + fitsPerChipWithLastPlanes[0]->XZ.error2At(quadCenter) );
				currentEntry.centerErrorPerChipFit.z = sqrt( fitsPerChipWithFirstPlanes[1]->YZ.error2At(quadCenter) + fitsPerChipWithLastPlanes[0]->YZ.error2At(quadCenter) );
			} else if(fitsPerChip[2] && fitsPerChip[3]) {
				currentEntry.centerDiffPerChipFit.x = fitsPerChipWithFirstPlanes[2]->xAt(quadCenter) - fitsPerChipWithLastPlanes[3]->xAt(quadCenter);
				currentEntry.centerDiffPerChipFit.z = fitsPerChipWithFirstPlanes[2]->yAt(quadCenter) - fitsPerChipWithLastPlanes[3]->yAt(quadCenter);
				currentEntry.centerErrorPerChipFit.x = sqrt( fitsPerChipWithFirstPlanes[2]->XZ.error2At(quadCenter) + fitsPerChipWithLastPlanes[3]->XZ.error2At(quadCenter) );
				currentEntry.centerErrorPerChipFit.z = sqrt( fitsPerChipWithFirstPlanes[2]->YZ.error2At(quadCenter) + fitsPerChipWithLastPlanes[3]->YZ.error2At(quadCenter) );
			} else {
				currentEntry.centerDiffPerChipFit.x = 0;
				currentEntry.centerDiffPerChipFit.z = 0;
				currentEntry.centerErrorPerChipFit.x = 0;
				currentEntry.centerErrorPerChipFit.z = 0;
			}

			currentEntry.nHitsPerChip=nHitsPerChip;
			currentEntry.nHitsPerChipValid=countHitsPerChip(quadHits,true);
			//currentEntry.quadHits=quadHits;
			currentEntry.matched=matched;

			currentEntry.telescopeTime=telescopeFitter.timestamp;
			currentEntry.triggerToA=quadFitter.reader.triggerToA;

			outputTree->Fill();
		}


		if(drawEvent && matched) {
			//		SimpleDetectorConfiguration setupForDrawing { 0,30 /*x*/, 0,42 /*y beam*/, -20,20/*z drift*/};
			auto setupForDrawing=simpleDetectorFromChipCorners(alignment.getAllChipCorners());
			setupForDrawing.minz=-1, setupForDrawing.maxz=11;
			//		SimpleDetectorConfiguration setupForDrawing { 10,40 /*x*/, 0,400 /*y beam*/, 0,40/*z drift*/};

//			std::cout<<"draw telescope hits\n";
//			telescopeFitter.drawEvent(telescopeHits,telescopeFits);
	//		HoughTransformer::drawClusters(telescopeHits, setupForDrawing);

			std::cout<<"draw "<< (matched?"matched":"") <<" quad hits (selected) "<<quadHits.size()<<" ("<<countTotalValidHits(quadHits)<<")\n";

			//3D
			const bool draw3D=true;
			if(draw3D) {
				HoughTransformer::drawCluster(quadHits,setupForDrawing);
				if(fittedTrack) fittedTrack->draw( setupForDrawing.ymin(), setupForDrawing.ymax() );
//				for (auto& f : telescopeFits)
//					f.draw( setupForDrawing.ymin(), setupForDrawing.ymax() );
//				for (auto& f : timepixFits)
//					f.draw( setupForDrawing.ymin(), setupForDrawing.ymax(), kTeal );
				drawQuadOutline(alignment, setupForDrawing.zmax() );
				gPad->GetPrimitive("eventDisplayLegend")->Draw(); //move legend to top again
				gPad->Update();
			}
			const bool draw2D=false;
			if(draw2D){
				//2D
				drawCluster2D(quadHits,setupForDrawing);
				alignment.drawChipEdges();
				for (auto& f : telescopeFits)
					f.XZ.draw( setupForDrawing.ymin(), setupForDrawing.ymax() );
				for (auto& f : timepixFits)
					f.XZ.draw( setupForDrawing.ymin(), setupForDrawing.ymax(), kTeal );
				gPad->Update();
			}

			gPad->Update();
			if(processDrawSignals()) break;

//			if(not pdfIsOpen) {
//				gPad->Print("eventDisplays.pdf(");
//				pdfIsOpen=true;
//			} else {
//				gPad->Print("eventDisplays.pdf");
//			}
		}

	}
	if(keepStatus) timepixStatusKeepers.writeBuffer([](StatusKeeper&s){s.reset();});

	std::cout<<"\ndone!\n";

	if(not drawEvent and outputTree) {
		static TCanvas* xCorrelation= new TCanvas();
		xCorrelation->cd();
		outputTree->Draw("XZ.intercept+XZ.slope*172:x", "", "", 1E4);
		static TCanvas* yCorrelation= new TCanvas();
		yCorrelation->cd();
		outputTree->Draw("YZ.intercept+YZ.slope*172:z", "fabs(XZ.intercept+XZ.slope*172-x)<3", "", 1E4);
	}

	if(pdfIsOpen) {
		gPad->Print("eventDisplays.pdf]");
	}
}

void TrackCombiner::printTriggers(int telescopeEntry, int tpcEntry) {
	int printEveryN=1;
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

		if(keepStatus) frameStatusHistogram->reset();
		nTelescopeTriggers+=std::max(0,telescopeFitter.triggerNumberEnd-previousTriggerNumberEnd);

		previousTriggerNumberBegin=telescopeFitter.triggerNumberBegin;
		previousTriggerNumberEnd=telescopeFitter.triggerNumberEnd;
		if( !telescopeFitter.getEntry(++telescopeEntry) ) return MatchResult::end;

		if(printMatching)		std::cout<<"increased telescopeEntry, first timepix match was: "<<timepixEntryFirstMatch<<std::endl;

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
