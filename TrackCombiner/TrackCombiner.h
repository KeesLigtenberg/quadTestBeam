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
#include "EntryBuffer.h"


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

//	const int triggerOffset=1000;//WARNING OUT OF SYNC!!!
	const int triggerOffset=0;
	int nTelescopeTriggers=0;
	int previousTriggerNumberBegin=0, previousTriggerNumberEnd=0;
	bool hadFirstMatch=false;
	int timepixEntryFirstMatch=0;

	enum class MatchResult { match, noMatch, end };
	TrackCombiner::MatchResult getAndMatchEntries(
			int& telescopeEntry,
			int& tpcStartEntry);




	//for making cut flow:
	struct StatusKeeper{
		StatusKeeper(std::string name) : statusHistogram( new TH1D( (name+"Status").c_str(), ("status of "+name).c_str(), 1,0,1) ) {};
		StatusKeeper(std::shared_ptr<TH1D> h) : statusHistogram(h) {};
		StatusKeeper(const TrackCombiner::StatusKeeper& o) : statusHistogram(o.statusHistogram) {};
		virtual ~StatusKeeper() {};
		int priority=0; std::string message="";
		void replace(int messagePriority, std::string newMessage) {
			if(messagePriority>priority) { message=newMessage; priority=messagePriority;}
		}
		void reset() {
			if(priority) statusHistogram->Fill(message.c_str(),1);
			priority=0;
		}
		void Write() { statusHistogram->LabelsDeflate(); statusHistogram->Write(); };
		std::shared_ptr<TH1D> statusHistogram;
	};
	const bool keepStatus=true;
	std::unique_ptr<StatusKeeper> frameStatusHistogram{}, triggerStatusHistogram{}, timepixStatusHistogram{};
	EntryBuffer<int, StatusKeeper> timepixStatusKeepers; //special buffered statusKeeper function

	//one function to fill all statushistograms at once.
	void replaceStatus(int priority, std::string message, int tpcEntryNumber) {
		if(not keepStatus) return;
		for(auto* s : {&frameStatusHistogram, &triggerStatusHistogram} ) (*s)->replace(priority, message);
		if(!timepixStatusKeepers.isInBuffer(tpcEntryNumber)) {
			timepixStatusKeepers.placeInBuffer(tpcEntryNumber, *timepixStatusHistogram);
		}
		timepixStatusKeepers.getFromBuffer(tpcEntryNumber).replace(priority, message);
	}


};

#pragma link C++ class std::vector<int>+;

#endif /* TRACKCOMBINER_TRACKCOMBINER_H_ */
