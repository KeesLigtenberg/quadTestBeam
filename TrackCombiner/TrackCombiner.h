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
	void openFile(std::string outputFileName);


	void printTriggers(int telescopeEntry, int tpcEntry);

	QuadTrackFitter quadFitter;
	TelescopeTrackFitter telescopeFitter;
	Alignment& alignment;

private:
	std::unique_ptr<TFile> outputFile{nullptr};
	std::unique_ptr<TTree> outputTree{nullptr};

	struct TreeEntry {
		std::vector<FitResult3D> telescopeFits;
		Vec3 meanQuadPosition;
		std::vector<int> nHitsPerChip;
		std::vector<PositionHit> quadHits;
		bool matched;
	} currentEntry;

	int triggerOffset=0;
	int nTelescopeTriggers=0;
	int previousTriggerNumberBegin=0, previousTriggerNumberEnd=0;
	bool hadFirstMatch=false;
	int timepixEntryFirstMatch=0;

	enum class MatchResult { match, noMatch, end };
	TrackCombiner::MatchResult getAndMatchEntries(
			int& telescopeEntry,
			int& tpcStartEntry);

};

#pragma link C++ class std::vector<int>+;

#endif /* TRACKCOMBINER_TRACKCOMBINER_H_ */
