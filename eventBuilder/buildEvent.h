/*
 * buildEvent.h
 *
 *  Created on: Oct 8, 2018
 *      Author: cligtenb
 */

#ifndef EVENTBUILDER_BUILDEVENT_H_
#define EVENTBUILDER_BUILDEVENT_H_


#include <string>
#include <array>
#include <vector>
#include <memory>
#include <algorithm>

#include "TTree.h"
#include "TH1.h"

#include "/user/cligtenb/rootmacros/getObjectFromFile.h"

#include "TriggerDecoder.h"
#include "Hit.h"

//run name from folder name
std::string getRunFromFolder(std::string runDir) {
	auto pos=runDir.rfind('/');
	if(pos!=std::string::npos) {//found
		return runDir.substr(pos);
	} else {
		return runDir;
	}
}

struct TreeReader {
	TTree* tree{};
	unsigned currentEntry{0};

	TreeReader() {};
	TreeReader(TTree* tree) : tree(tree) {};

	virtual ~TreeReader() {};

	virtual void getFirstEntry() {
		currentEntry=0;
		tree->GetEntry(currentEntry);
	}

	virtual void getNextEntry() {
		tree->GetEntry(++currentEntry);
	}
	bool reachedEnd() {
		return currentEntry+1>=tree->GetEntries();
	}
};


struct ChipTreeReader : TreeReader {
	unsigned long long toa{};
	unsigned char col{},row{};
	unsigned short tot{};

	ChipTreeReader( TTree* tree ) : TreeReader(tree) {
		tree->SetBranchAddress("toa", &toa);
		tree->SetBranchAddress("col", &col);
		tree->SetBranchAddress("row", &row);
		tree->SetBranchAddress("tot", &tot);
		tree->GetEntry(currentEntry);
	}
	void discardEntriesBeforeT(unsigned long long t) {
		while(not reachedEnd() and toa < t) getNextEntry();
	}
};

struct TriggerTreeReader : TreeReader {
	unsigned long long toa{};
	TriggerTreeReader( TTree* tree ) : TreeReader(tree) {
		tree->SetBranchAddress("timestamp", &toa);
		tree->GetEntry(currentEntry);
	}
	virtual void getNextEntry() {
		auto oldToa=toa;
		tree->GetEntry(++currentEntry);
		while(oldToa==toa) {
			static int warnedDoubleTrigger=0;
			warnedDoubleTrigger++;
			if(warnedDoubleTrigger<20) std::cout<<"Warning: double entry in trigger tree at "<<currentEntry<<"\n";
			if(warnedDoubleTrigger==20) std::cout<<"Suppressing further warnings on double entries in trigger tree\n";
			tree->GetEntry(++currentEntry);
		}
	}

	struct Trigger {
		unsigned long long toa;
		unsigned number;
	};
	std::vector<int> getNextTriggerWord( unsigned long long firstToa) {
		unsigned timeDifference=0;
		std::vector<int> triggerword;
		//read all bits in this triggernumber
		while(not reachedEnd() ){
				getNextEntry();
//			std::cout<<currentEntry<<" "<<toa<<" ";
			timeDifference=toa-firstToa;
			if( timeDifference<newTriggerThreshold) {
				triggerword.push_back(timeDifference);
//				std::cout<<"timeDifference: "<<timeDifference<<"\n";
			}
			else break;
		}
		return triggerword;
	}
	Trigger getNextTrigger() {
		Trigger trigger;
		trigger.toa=toa; //first edge is trigger time
		auto triggerword=getNextTriggerWord(trigger.toa);
		trigger.number = decoder.getNextTriggerFromTimes(triggerword);
		return trigger;
	}

	void updateOffsetAndBinWidth(int sampleSize=20) {
		getFirstEntry();
//		std::cout<<"got first entry\n";
		double firstBinOffsetSum=0, timeBinWidthSum=0;
		int firstBinOffsetEntries=0, timeBinWidthEntries=0;
		for(int i=1; i<=sampleSize; i++) {
			auto firstToa=toa;
			auto triggerword=getNextTriggerWord(firstToa);

			if(!triggerword.size() and i==1) {
				std::cout<<"first trigger is 0\n";
				i=0; continue;
			}

			if(i%2) { //triggernumber is odd?'
//				std::cout<<"trigger number "<<i<<" trying to find offset\n";
				firstBinOffsetSum+=triggerword.front();
				++firstBinOffsetEntries;
			}
			if(triggerword.size()>1) {
				std::vector<int> triggerDiffs;
				 for(unsigned i=0; i+1<triggerword.size(); i++) {
					triggerDiffs.push_back(triggerword.at(i+1)-triggerword.at(i) );
//					std::cout<<"diff: "<<triggerDiffs.back()<<"\n";
				}
				auto min=std::min_element(triggerDiffs.begin(), triggerDiffs.end());
				int usedEdge=2;
				if(*min>(usedEdge-0.5)*decoder.binTimeWidth and *min<(usedEdge+0.5)*decoder.binTimeWidth) {
					timeBinWidthSum+=*min/usedEdge;
					++timeBinWidthEntries;
				}
			}
		}

		std::cout<<"updated binwidth "<<decoder.binTimeWidth<<" and first bin offset "<<decoder.firstBinOffset;
		decoder.binTimeWidth=timeBinWidthSum/timeBinWidthEntries;
		decoder.firstBinOffset=firstBinOffsetSum/firstBinOffsetEntries;
		std::cout<<"\nto binwidth "<<decoder.binTimeWidth<<" and first bin offset "<<decoder.firstBinOffset<<"\n\n";

		getFirstEntry();
	}

	const double newTriggerThreshold=5E3/25*4096; //time between two seperate triggers
	const int nbits=15;
	TriggerDecoder decoder{nbits, 200./25*4096/*timeBinWidth*/, 564/25.*4096 /*first bin offset*/};

};



struct CombinedTreeWriter {
	TTree tree;
	std::vector<Hit> chips[4];
	long long triggerToA;
	unsigned triggerNumber;

	CombinedTreeWriter() {
		tree.Branch( "triggerToA", &triggerToA );
		tree.Branch( "triggerNumber", &triggerNumber );
		for(int i=0; i<4; i++)
			tree.Branch( ("chip"+std::to_string(i)).c_str() , "std::vector<Hit>", chips+i);
	}
	void fill() {
			tree.Fill();
			for(auto& c : chips) c.clear();
	}
};

void convertToTree(std::string inputFileName, std::string outputFileName) {
	//Read trees
	TFile inputFile(inputFileName.c_str(), "READ");
	std::cout<<"reading trees\n";
	std::vector<std::unique_ptr<ChipTreeReader> > chips;
	for(int i=0; i<4; i++) {
		std::cout<<"  chip "<<i<<"\n";
		auto chipTree=getObjectFromFile<TTree>("hits_chip"+std::to_string(i), &inputFile);
		chips.emplace_back( new ChipTreeReader(chipTree) );
	}
	std::cout<<"  trigger\n";
	auto triggerTree=getObjectFromFile<TTree>("triggers", &inputFile);
	TriggerTreeReader triggerReader(triggerTree);

	//Output tree
	std::cout<<"creating output file\n";
	TFile outputFile((outputFileName).c_str(), "RECREATE"); //use absolute path
	CombinedTreeWriter outputTrees;

	std::cout<<"updating offset and bin width..\n";
//	triggerReader.updateOffsetAndBinWidth(200);

	while( not triggerReader.reachedEnd()) {
		TriggerTreeReader::Trigger trigger;
		try{
			trigger=triggerReader.getNextTrigger();
		} catch (TriggerDecoder::TriggerDecoderException& e) {
			if(e.description=="wrong zero!") { std::cout<<e.description<<"\n"; continue; }
			else throw e;
		}
		if(!(trigger.number%10000)) {
			std::cout<<trigger.number<<" at time "<<trigger.toa<<"\n";
		}
//		if(trigger.number>50000) break;

		//set triggerReader
		outputTrees.triggerToA=trigger.toa;
		outputTrees.triggerNumber=trigger.number;

		//find all hits within a range of the triggerReader
		const unsigned maxTimeBeforeTrigger=500/25*4096; //500 ns
		const unsigned maxTimeAfterTrigger=500/25*4096; //500 ns
		for(unsigned i=0; i<chips.size(); i++) {
			auto& c = chips[i];
			c->discardEntriesBeforeT(trigger.toa-maxTimeBeforeTrigger);
			while(c->toa < trigger.toa+maxTimeAfterTrigger and not c->reachedEnd()) {
				outputTrees.chips[i].emplace_back(c->row, c->col, c->tot, c->toa-trigger.toa );
				c->getNextEntry();
			}
		}//for chips

		//fill tree after each triggerReader
		outputTrees.fill();
	}

	outputTrees.tree.Write("data");
	outputFile.Write();
}


void buildEvent() {





}




#endif /* EVENTBUILDER_BUILDEVENT_H_ */
