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

	bool drawEvent=false;

private:
	std::unique_ptr<TFile> outputFile{nullptr};
	std::unique_ptr<TTree> outputTree{nullptr};
	std::array< std::unique_ptr<ChipHistogrammer>, nChips> hists{{}};
	std::unique_ptr<ChipHistogrammer> quadHist{nullptr};
	std::unique_ptr<TH1D> selectedHitAverageToTrackx, selectedHitAverageToTrackz, fractionInTrack, smallestShift, averageShift;
	std::unique_ptr<TH1D> averageToTrackx;


	struct TreeEntry {
		std::vector<FitResult3D> telescopeFits;
		std::vector<FitResult3D> timepixFits;
		Vec3 meanQuadPosition, meanQuadDiff, meanQuadError;
		std::vector<Vec3> meanPositionPerChip, meanDiffPerChip, meanErrorPerChip;
		std::vector<std::vector<Vec3>> meanDiffPerChipPerFitFirst, meanDiffPerChipPerFitLast;
		std::vector<std::vector<Vec3>> meanErrorPerChipPerFitFirst, meanErrorPerChipPerFitLast;
		Vec3 centerDiffPerChipFit, centerErrorPerChipFit;
		std::vector<int> nHitsPerChip;
		std::vector<int> nHitsPerChipValid;
		std::vector<PositionHit> quadHits;
		bool matched;
		long long triggerToA;
		unsigned telescopeTime;
	} currentEntry;

//	const int triggerOffset=100;//WARNING OUT OF SYNC!!!
	const int triggerOffset=0;
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
