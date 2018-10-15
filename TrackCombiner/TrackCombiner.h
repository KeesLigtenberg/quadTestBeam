/*
 * TrackCombiner.h
 *
 *  Created on: Oct 15, 2018
 *      Author: cligtenb
 */

#ifndef TRACKCOMBINER_TRACKCOMBINER_H_
#define TRACKCOMBINER_TRACKCOMBINER_H_

#include "../TelescopeTrackFitter/TelescopeTrackFitter.cpp"
#include "../TrackFitter/QuadTrackFitter.cpp"

class TrackCombiner {
public:
	TrackCombiner(std::string quadFile, std::string telescopeFile, Alignment& alignment);
	virtual ~TrackCombiner();

	void Process();

private:

	void printTriggers(int telescopeEntry, int tpcEntry) const;

	TelescopeTrackFitter telescopeFitter;
	QuadTrackFitter quadFitter;
	Alignment& alignment;


	int triggerOffset=-1;
	int nTelescopeTriggers=0;
	int previousTriggerNumberBegin=0, previousTriggerNumberEnd=0;
	bool hadFirstMatch=false;
	int timepixEntryFirstMatch=0;

	enum class MatchResult { match, noMatch, end };
	TrackCombiner::MatchResult TrackCombiner::getAndMatchEntries(
			int& telescopeEntry,
			int& tpcStartEntry);

};

#endif /* TRACKCOMBINER_TRACKCOMBINER_H_ */
