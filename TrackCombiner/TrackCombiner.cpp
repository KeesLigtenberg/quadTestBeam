/*
 * TrackCombiner.cpp
 *
 *  Created on: Oct 15, 2018
 *      Author: cligtenb
 */

#include "TrackCombiner.h"

#include <memory>


namespace {

	PositionHit& calculateResidualTimepix(PositionHit& h, const FitResult3D& telescopeFit ) {
		h.residual.x=h.position.x - telescopeFit.xAt(h.position.y);
		h.residual.z=h.position.z - telescopeFit.yAt(h.position.y);
		h.residual.y=0;
		return h;
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
	outputTree=std::unique_ptr<TTree>(new TTree("fitResults", "fitResults"));

	outputTree->Branch("telescopeFits", "std::vector<FitResult3D>", &currentEntry.telescopeFits);
	outputTree->Branch("meanQuadPosition", "Vec3", &currentEntry.meanQuadPosition);
	outputTree->Branch("nHitsPerChip", "std::vector<int>", &currentEntry.nHitsPerChip);
	outputTree->Branch("quadHits", "std::vector<PositionHit>", &currentEntry.quadHits);
	outputTree->Branch("matched", &currentEntry.matched);

	for(int i=0; i<nChips; i++) {
		hists[i]=std::unique_ptr<ChipHistogrammer>( new ChipHistogrammer("chip"+std::to_string(i), alignment) );
	}
	quadHist=std::unique_ptr<ChipHistogrammer>( new ChipHistogrammer("quad", alignment) );

}

void TrackCombiner::Process() {
	//match entries
	nTelescopeTriggers=0;
	telescopeFitter.getEntry(previousTriggerNumberBegin);

	quadFitter.reader.tree->GetEntry(0);


	for(int telescopeEntryNumber=0,tpcEntryNumber=0; //5000000, 2308829
			telescopeEntryNumber<1E4 //telescopeFitter.nEvents//1000000
			;) {
//		triggerStatusHistogram.reset();

		// Get Entry and match trigger Numbers
		auto matchStatus=getAndMatchEntries(telescopeEntryNumber,tpcEntryNumber);
		printTriggers(telescopeEntryNumber,tpcEntryNumber);
		if( matchStatus == MatchResult::end ) {
//			std::cout<<"match is end\n\n";
			break;}
		else if( matchStatus == MatchResult::noMatch) {
//			std::cout<<" no match\n\n";
			continue; }
//		else { std::cout<<"match!\n\n"; }


		auto telescopeHits=telescopeFitter.getSpaceHits();
		int nTelescopeHits=0;
		for(const auto& plane : telescopeHits) nTelescopeHits+=plane.size();
//		std::cout<<"there are "<<nTelescopeHits<<" telescope hits\n";
		if(!nTelescopeHits) continue;

 		auto quadHits=quadFitter.getSpaceHits(alignment);
 		if(!quadHits.size()) continue;
 		auto nHitsPerChip= getHitsPerChip(quadHits);

 		const int minHitsPerChip=10;
 		if( (nHitsPerChip[0] <minHitsPerChip or nHitsPerChip[1] <minHitsPerChip) and (nHitsPerChip[2]<minHitsPerChip or nHitsPerChip[3]<minHitsPerChip ) ) continue;


 		std::vector<FitResult3D> telescopeFits=telescopeFitter.getFits(telescopeHits);


 		auto avPosition=getAveragePosition(quadHits);
 		bool matched=false;
 		for(auto&f : telescopeFits) {
 			if( std::fabs(avPosition.x()-f.xAt(193.25))<2 and std::fabs(avPosition.z()-f.yAt(193.25))<2 ) {
 				matched=true;
 				for(auto& h : quadHits) {
 					h=calculateResidualTimepix(h, f);
 				}
 			}
 		}
// 		if(!matched) continue;


 		if(matched) {
 			for(auto& h : quadHits) {
 				hists[h.chip]->fillTimewalk(h, alignment.timeWalk);
 				quadHist->fillTimewalk(h, alignment.timeWalk);

 				hists[h.chip]->fillHit(h);
 				hists[h.chip]->fillRotation(h, alignment.chips[h.chip].COM);
 				quadHist->fillHit(h);
 				quadHist->fillRotation(h, alignment.quad.COM);
 			}
 		}


		for(auto& h : hists) {
			h->fillEvent();
		}
		quadHist->fillEvent();

		const bool drawEvent=true;
		if(drawEvent) {
			//		SimpleDetectorConfiguration setupForDrawing { 0,30 /*x*/, 0,42 /*y beam*/, -20,20/*z drift*/};
			auto setupForDrawing=simpleDetectorFromChipCorners(alignment.getAllChipCorners());
			setupForDrawing.minz=-10, setupForDrawing.maxz=30;
			//		SimpleDetectorConfiguration setupForDrawing { 10,40 /*x*/, 0,400 /*y beam*/, 0,40/*z drift*/};

//			std::cout<<"draw telescope hits\n";
//			telescopeFitter.drawEvent(telescopeHits,telescopeFits);
	//		HoughTransformer::drawClusters(telescopeHits, setupForDrawing);

//			std::cout<<"draw quad hits "<<quadHits.size()<<"\n";
//			HoughTransformer::drawCluster(quadHits,setupForDrawing);
//			for (auto& f : telescopeFits)
//				f.draw( setupForDrawing.ymin(), setupForDrawing.ymax() );
//			drawQuadOutline(alignment, setupForDrawing.zmin() );

			drawCluster2D(quadHits,setupForDrawing);
			alignment.drawChipEdges();
			for (auto& f : telescopeFits)
				f.XZ.draw( setupForDrawing.ymin(), setupForDrawing.ymax() );


			gPad->Update();
			if(processDrawSignals()) break;
		}

		if(outputTree) {
			currentEntry.telescopeFits=telescopeFits;
			currentEntry.meanQuadPosition=getAveragePosition(quadHits);
			currentEntry.nHitsPerChip=nHitsPerChip;
			currentEntry.quadHits=quadHits;
			currentEntry.matched=matched;
			outputTree->Fill();
		}

	}
	outputTree->Draw("XZ.intercept+XZ.slope*172:x");
	new TCanvas();
	outputTree->Draw("YZ.intercept+YZ.slope*172:z", "fabs(XZ.intercept+XZ.slope*172-x)<3");
}

void TrackCombiner::printTriggers(int telescopeEntry, int tpcEntry) {
	const int printEveryN=1000;
	if( !(telescopeEntry%printEveryN) ) {
		cout<<"entry: "<<telescopeEntry<<"/"<<telescopeFitter.nEvents<<" ";
		cout<<"triggers: "<<telescopeFitter.triggerNumberBegin<<"-"<<telescopeFitter.triggerNumberEnd;
		cout<<" timepix triggerNumber: "<<quadFitter.triggerNumber()<<"="<<(quadFitter.triggerNumber()+triggerOffset) % 32768<<" in entry "<<tpcEntry<<endl;
	}
}

TrackCombiner::MatchResult TrackCombiner::getAndMatchEntries(
		int& telescopeEntry,
		int& tpcStartEntry) {

//	std::cout<<"starting getAndMatchEntries( "<<telescopeEntry<<", "<<tpcStartEntry<<")\n";

	//use modulus offset instead of doing actual modulus, to continue counting after 32768
	long long modulusOffset=32768*int((quadFitter.triggerNumber()+triggerOffset)/32768);

	//special case: triggerNumberBegin==triggerNumberEnd==0
	while(telescopeFitter.triggerNumberBegin==0 and telescopeFitter.triggerNumberEnd==0) {
		std::cout<<"special case telescope trigger number begin==end==0\n";
		if( !telescopeFitter.getEntry(++telescopeEntry) ) return MatchResult::end;
	}

	//if telescope triggerNumberBegin decreased, and tpc did not already pass this boundary
	//then we must get entries until the tpc triggernumber also passes the 32768 boundary ( is larger or equal to new trigger number)
	if(telescopeFitter.triggerNumberBegin<previousTriggerNumberBegin
			and not (quadFitter.triggerNumber()+triggerOffset-modulusOffset<500) ) {
//		std::cout<<"telescope number decreased and tpc did not yet, get entries until tpc also passes boundary\n";
		do {
//			printTriggers(telescopeEntry, tpcStartEntry);
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
//		std::cout<<"getting tpc, until tpc trigger entry is larger or equal to begin\n";
		if( not quadFitter.getEntry(tpcStartEntry++) ) {
			std::cout<<"could not get entry "<<tpcStartEntry-1<<"\n";
			return MatchResult::end;
		}
	} while( (quadFitter.triggerNumber()+triggerOffset-modulusOffset) < telescopeFitter.triggerNumberBegin
			and telescopeFitter.triggerNumberBegin - (quadFitter.triggerNumber()+triggerOffset-modulusOffset)<500 );

	//if also larger than end: reached the end of this telescope frame, continue with next telescope frame;
	if( (quadFitter.triggerNumber()+triggerOffset-modulusOffset) > telescopeFitter.triggerNumberEnd) {
//		std::cout<<"continuing to next telescope frame\n";
//		triggerStatusHistogram.Fill("Trigger numbers do not match", 1);

//		frameStatusHistogram.reset();
		nTelescopeTriggers+=max(0,telescopeFitter.triggerNumberEnd-previousTriggerNumberEnd);

		previousTriggerNumberBegin=telescopeFitter.triggerNumberBegin;
		previousTriggerNumberEnd=telescopeFitter.triggerNumberEnd;
		if( !telescopeFitter.getEntry(++telescopeEntry) ) return MatchResult::end;

//		cout<<"increased telescopeEntry, first time pix match was: "<<timepixEntryFirstMatch<<endl;

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
