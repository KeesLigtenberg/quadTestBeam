/*
 * TrackCombiner.cpp
 *
 *  Created on: Oct 15, 2018
 *      Author: cligtenb
 */

#include "TrackCombiner.h"

TrackCombiner::TrackCombiner(std::string quadFile, std::string telescopeFile, Alignment& alignment) :
	quadFitter(quadFile),
	telescopeFitter(telescopeFile, mimosa),
	alignment(alignment)
{
	telescopeFitter.setAlignment(alignment);
	telescopeFitter.makeMask(1E4);
}

TrackCombiner::~TrackCombiner() {
}

void TrackCombiner::Process() {

	//match entries
	nTelescopeTriggers=0;
	telescopeFitter.getEntry(previousTriggerNumberBegin);
	for(int telescopeEntryNumber=0,tpcEntryNumber=0; //5000000, 2308829
			telescopeEntryNumber<2E6//telescopeFitter.nEvents//1000000
			;) {
//		triggerStatusHistogram.reset();

		// Get Entry and match trigger Numbers
		auto matchStatus=getAndMatchEntries(telescopeEntryNumber,tpcEntryNumber);
		printTriggers(telescopeEntryNumber,tpcEntryNumber);
//		if( cin.get()=='q') break;
//		if(!(telescopeEntryNumber%100000)) cout<<"entry "<<telescopeEntryNumber<<endl;
		if( matchStatus == MatchResult::end ) break;
		else if( matchStatus == MatchResult::noMatch) continue;


		auto telescopeHits=telescopeFitter.getSpaceHits();
		auto quadHits=quadFitter.getSpaceHits(alignment);


		SimpleDetectorConfiguration setupForDrawing { 10,40 /*x*/, 0,42 /*y beam*/, 0,40/*z drift*/};
		HoughTransformer::drawCluster(quadHits,setupForDrawing);
		if(processDrawSignals()) break;

	}
}

void TrackCombiner::printTriggers(int telescopeEntry, int tpcEntry) const {
	const int printEveryN=10000;
	if( !(telescopeEntry%printEveryN) ) {
		cout<<"entry: "<<telescopeEntry<<"/"<<telescopeFitter.nEvents<<" ";
		cout<<"triggers: "<<telescopeFitter.triggerNumberBegin<<"-"<<telescopeFitter.triggerNumberEnd;
		cout<<" timepix triggerNumber: "<<quadFitter.triggerNumber()<<"="<<(quadFitter.triggerNumber()+triggerOffset) % 32768<<" in entry "<<tpcEntry<<endl;
	}
}

TrackCombiner::MatchResult TrackCombiner::getAndMatchEntries(
		int& telescopeEntry,
		int& tpcStartEntry) {

	//use modulus offset instead of doing actual modulus, to continue counting after 32768
	long long modulusOffset=32768*int((quadFitter.triggerNumber()+triggerOffset)/32768);

	//if telescope triggerNumberBegin decreased, and tpc did not already pass this boundary
	//then we must get entries until the tpc triggernumber also passes the 32768 boundary ( is larger or equal to new trigger number)
	if(telescopeFitter.triggerNumberBegin<previousTriggerNumberBegin
			and not (quadFitter.triggerNumber()+triggerOffset-modulusOffset<500) ) {
		do {
			if( !quadFitter.getEntry(tpcStartEntry++) ) return MatchResult::end;
		} while ( quadFitter.triggerNumber()+triggerOffset - modulusOffset - 32768 < telescopeFitter.triggerNumberBegin );
		modulusOffset=32768*int((quadFitter.triggerNumber()+triggerOffset)/32768);
		//stop timepix from going back over this boundary
		timepixEntryFirstMatch=tpcStartEntry;
		hadFirstMatch=true;
	}

	//get next entry until tpc trigger number is larger than or equal to begin
	//or until the tpc trigger number decreases;
	do {
		if( !quadFitter.getEntry(tpcStartEntry++) ) return MatchResult::end;
	} while( (quadFitter.triggerNumber()+triggerOffset-modulusOffset) < telescopeFitter.triggerNumberBegin
			and telescopeFitter.triggerNumberBegin - (quadFitter.triggerNumber()+triggerOffset-modulusOffset)<500 );

	//if also larger than end: reached the end of this telescope frame, continue with next telescope frame;
	if( (quadFitter.triggerNumber()+triggerOffset-modulusOffset) > telescopeFitter.triggerNumberEnd) {
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
