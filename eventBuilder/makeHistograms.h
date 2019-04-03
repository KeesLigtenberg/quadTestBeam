/*
 * makeHistograms.h
 *
 *  Created on: Oct 31, 2018
 *      Author: cligtenb
 */

#ifndef EVENTBUILDER_MAKEHISTOGRAMS_H_
#define EVENTBUILDER_MAKEHISTOGRAMS_H_

#include "buildEvent.h"

#include <deque>
#include <list>

#include "TH1.h"
#include "TH2.h"
#include "TProfile2D.h"
#include "TCanvas.h"
#include "TStyle.h"

void makeHistograms(std::string inputFileName, std::string outputFileName) {
	//Read trees
	TFile inputFile(inputFileName.c_str(), "READ");
	std::cout<<"reading trees\n";
	const int nChips=4;
	std::vector<std::unique_ptr<ChipTreeReader> > chips;
	for(int i=0; i<nChips; i++) {
		std::cout<<"  chip "<<i<<"\n";
		auto chipTree=getObjectFromFile<TTree>("hits_chip"+std::to_string(i), &inputFile);
		chips.emplace_back( new ChipTreeReader(chipTree) );
		chips[i]->getFirstEntry();
	}
	std::cout<<"  trigger\n";
	auto triggerTree=getObjectFromFile<TTree>("triggers", &inputFile);
	TriggerTreeReader triggerReader(triggerTree);

	//Output tree
	std::cout<<"creating output file\n";
	TFile outputFile((outputFileName).c_str(), "RECREATE"); //use absolute path
	CombinedTreeWriter outputTrees;

	//histograms
	TH2D* pixelHitMap[nChips]{};
	for(int i=0; i<nChips; i++) pixelHitMap[i]=new TH2D(("chip"+std::to_string(i)).c_str(), "pixel hitmap", 256,0,256,256,0,256);

	auto canv=new TCanvas("hitmap", "hitmap", 1000, 1000);
	canv->Divide(2,2);
	gStyle->SetPalette(1);

	long frame=0;

	std::deque<TriggerTreeReader::Trigger> triggers;
	int triggerWindowUs=8E3;
	int stepSizeUs=5; //409; //5
	unsigned triggerWindow=triggerWindowUs*1E3/25*4096; //to stepsize units
	TH1D triggerTime("triggerTime", "Trigger time offset;t - trigger time [us];Entries", triggerWindowUs/stepSizeUs, -triggerWindowUs/2,triggerWindowUs/2);
	std::vector<std::unique_ptr<TH1D>> triggerTimeOffsetChip;
	for(int i=0; i<nChips; i++) {
		triggerTimeOffsetChip.emplace_back( new TH1D(("triggerTimeOffsetChip"+std::to_string(i)).c_str(), "Trigger time offset;t - trigger time [us];Entries", triggerWindowUs/stepSizeUs, -triggerWindowUs/2,triggerWindowUs/2 ) );
	}
	const int nBits=64;
	TProfile2D triggerBitsByTime("triggerBitsByTime", "Trigger time offset;t - trigger time [us];trigger bits", triggerWindowUs/stepSizeUs, -triggerWindowUs/2,triggerWindowUs/2, nBits, 0.5,nBits+.5);
	TProfile2D hitTimeBitsByTime("hitTimeBitsByTime", "Hit time time offset;t - trigger time [us];trigger bits", triggerWindowUs/stepSizeUs, -triggerWindowUs/2,triggerWindowUs/2, nBits, 0.5,nBits+.5);
	TH2D triggerIndexByTime("triggerIndexByTime", "Hit time time offset;t - trigger time [us];trigger index", triggerWindowUs/stepSizeUs, -triggerWindowUs/2,triggerWindowUs/2, 100, 0,100);
	TH1D timeBetweenTriggers("timeBetweenTriggers", "Time to previous trigger; Time to previous trigger [us]; Entries", triggerWindowUs/stepSizeUs, 0, 1000);
	TH1D timeBetweenTracks("timeBetweenTracks", "Time to previous track (hits>N); Time to previous track [us]; Entries", triggerWindowUs/stepSizeUs, 0, 1000);
	TH1D nHitsTotal("nHitsTotal", "Total number of hits", 1000,0,1000);
	std::vector<std::unique_ptr<TH1D> > nHitsPerChip;
	for(int i=0; i<4; i++) {
		nHitsPerChip.emplace_back( new TH1D( ("nHitsPerChip"+std::to_string(i)).c_str(), ("Number of hits on chip "+std::to_string(i)).c_str(), 1000, 0, 1000)  ) ;
	}

//	std::cout<<"updating offset and bin width..\n";
//	triggerReader.updateOffsetAndBinWidth(200);

	auto trigger=triggerReader.getNextTriggerForced();

	//move all chips to first trigger entry
	for(unsigned i=0; i<chips.size(); i++) {
		chips[i]->discardEntriesBeforeT(trigger.toa);
	}

	auto startToA=(*std::min_element(chips.begin(),chips.end(), [](std::unique_ptr<ChipTreeReader>& tr1, std::unique_ptr<ChipTreeReader>& tr2){return tr1->toa<tr2->toa;}))->toa;
	unsigned stepSize=stepSizeUs*1E3/25*4096; //to stepsize units
//	startToA+=1E9/25*4096;
	double previousTrackt=0;
//	std::cout<<"start is "<<startToA<<"\n";
	for(auto t=startToA; not chips[0]->reachedEnd() and not triggerReader.reachedEnd(); t+=stepSize) {
		if(!(frame%10000)) {
			std::cout<<"frame "<<frame<<" saved "<<outputTrees.triggerNumber<<"\r"<<std::flush;
			if(!(frame%1000000)) {  std::cout<<"\nforcing autosave..\n"; outputTrees.tree.AutoSave();}
		}
//		if( frame > 10E6  ) break;
		if(outputTrees.triggerNumber>1E5) break;

		//read triggers
		while(triggers.size() and triggers.front().toa< t -triggerWindow/2) triggers.pop_front();
		while (triggers.empty() or trigger.toa<t+triggerWindow/2) {
			auto previousTrigger=trigger;
			trigger=triggerReader.getNextTriggerForced();
			timeBetweenTriggers.Fill( (trigger.toa-previousTrigger.toa)*25E-3/4096 );
			triggers.push_back(trigger);
		}
//		std::cout<<"t="<<t<<" front trigger toa="<<triggers.front().toa<<" added trigger toa="<<trigger.toa<<"\n";

		std::array<int,nChips> nHits{{}};
		outputTrees.triggerToA=t;
		for(unsigned i=0; i<chips.size(); i++) {
			auto& c = chips[i];
			c->discardEntriesBeforeT(t);
			while(c->toa < t+stepSize and not c->reachedEnd()) {
				pixelHitMap[i]->Fill(c->col, c->row);
				outputTrees.chips[i].emplace_back(c->row, c->col, c->tot, c->toa-t );
				c->getNextEntry();
				nHits[i]++;
			}
		}//for chips
		frame++;

//		if(true) {
		if(nHits[0]+nHits[1]>20 or nHits[2]+nHits[3]>20) {
//		if(nHits[0]>40) {
//			if(nHits[0]+nHits[1]>70) {
//				outputTrees.chips[2].clear();
//				outputTrees.chips[3].clear();
//			}
//			if(nHits[2]+nHits[3]>70) {
//				outputTrees.chips[0].clear();
//				outputTrees.chips[1].clear();
//			}

			//TEMPORARELY DISABLE TREE WRITING!!!
			//outputTrees.fill();
			outputTrees.triggerNumber++;

//			std::cout<<triggers.size()<<"\n";

			unsigned j=0;
			for(auto& iTrigger : triggers) {
//				std::cout<<(t-iTrigger.toa)*0.025/4096<<"\n";
				triggerTime.Fill( int(t-iTrigger.toa)*0.025/4096 ); //in us
				triggerIndexByTime.Fill( int(t-iTrigger.toa)*0.025/4096 , j++);


				auto triggerBits(iTrigger.toa);
				auto hitBits(t);
				for(int i=1; i<=nBits; i++) {
					hitTimeBitsByTime.Fill(int(t-iTrigger.toa)*0.025/4096, i, hitBits&0x1);
					triggerBitsByTime.Fill(int(t-iTrigger.toa)*0.025/4096, i, triggerBits&0x1);
					triggerBits>>=1;
					hitBits>>=1;
				}

//				if( fabs(t-iTrigger.toa) < stepSize ) { //remove trigger if a match at 0
//					iTrigger.toa=0;
//				}

			}
			timeBetweenTracks.Fill( int(t-previousTrackt)*.025/4096 );
			previousTrackt=t;

			if(std::cin.get()=='q') break;

			nHitsTotal.Fill(nHits[0]+nHits[1]+nHits[2]+nHits[3]);
			for(int i=0; i<4; i++) {
				if(nHits[i]) nHitsPerChip[i]->Fill(nHits[i]);
			}

		}
		for(int i=0; i<4; i++) {
			if(nHits[i]>20) {
				for(auto& iTrigger : triggers) {
	//				std::cout<<(t-iTrigger.toa)*0.025/4096<<"\n";
					triggerTimeOffsetChip[i]->Fill( int(t-iTrigger.toa)*0.025/4096 ); //in us
				}
			}
		}

		outputTrees.clear();

//		if(nHits>200){
//			for(int i=0; i<4;i++) {
//				canv->cd(i+1);
//				pixelHitMap[i]->Draw("colz");
//			}
//			canv->Update();
//			auto s=std::cin.get();
//			if(s=='q') break;
//			for(int i=0; i<4; i++) pixelHitMap[i]->Reset();
//		}
//		std::cout<<"\n";

	}
	std::cout<<"\n";

	std::cout<<"found "<<outputTrees.triggerNumber<<" saveable frames in "<<frame<<" frames\n";

	for(int i=0; i<nChips;i++) {
		canv->cd(i+1);
		pixelHitMap[i]->Draw("colz");
	}
	canv->Update();

	outputTrees.tree.Write();
	outputFile.Write();

}



#endif /* EVENTBUILDER_MAKEHISTOGRAMS_H_ */
