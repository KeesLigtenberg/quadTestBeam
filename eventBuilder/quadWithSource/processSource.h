//root macro

#include <string>
#include <iostream>
#include <cmath>
#include <algorithm>

#include "TPad.h"
#include "TCanvas.h"
#include "TProfile2D.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TRandom.h"

#include "../../TrackFitter/PositionHit.h"

#include "CrossTalkCalculator.h"


#pragma link C++ class std::vector<double>+;
#pragma link C++ class std::vector<std::pair<double, double> >+;

struct QuadTreeReader {
	QuadTreeReader(std::string treeName, std::string fileName) :
				file(new TFile(fileName.c_str(), "READ")),
				tree(nullptr) {
		tree= dynamic_cast<TTree*>(file->Get(treeName.c_str()));
		if(!tree) {
			std::cerr<<"could not get tree "<<treeName<<" from file "<<file->GetName()<<std::endl;
			throw 1;
		}

		tree->SetBranchAddress("triggerToA", &triggerToA);
		tree->SetBranchAddress("triggerNumber", &triggerNumber);
		for(int i=0; i<4; i++) {
			tree->SetBranchAddress( ("chip"+std::to_string(i)).c_str(), chip+i );
		}
	}

	Long64_t entryNumber{-1};
	bool getNext() {
		if(!(entryNumber%1000)) std::cout<<"entry "<<entryNumber<<"/"<<tree->GetEntries()<<"\n";
		if(++entryNumber>=tree->GetEntries()) return false;
		tree->GetEntry(entryNumber);
		return true;
	}

	std::unique_ptr<TFile> file;
	TTree* tree;

	Long64_t triggerToA{};
  unsigned triggerNumber{};
	std::vector<Hit>* chip[4]{};
};

struct SourceClusterWriter {
	SourceClusterWriter(std::string fileName) :
		file(new TFile(fileName.c_str(), "RECREATE")),
		tree(new TTree("clusters", "tree with clusters"))
	{
		tree->Branch("quadHits", "std::vector<PositionHit>", &hits);
		tree->Branch("nHits", &nHits);
		tree->Branch("clusterChip", &clusterChip);
		tree->Branch("meanCol", &meanCol);
		tree->Branch("meanRow", &meanRow);
		tree->Branch("meanToT", &meanToT);
		tree->Branch("meanToA", &meanToA);
		tree->Branch("stdDevCol", &stdDevCol);
		tree->Branch("stdDevRow", &stdDevRow);
		tree->Branch("frameToA", &frameToA);

		tree->Branch("nNeighbours",&nNeighbours);
		tree->Branch("isolatedToTs", "std::vector<double>", &isolatedToTs);
		tree->Branch("isolatedPairToTs", "std::vector<std::pair<double,double>>", &isolatedPairToTs);
		tree->Branch("controlPairToTs", "std::vector<std::pair<double,double>>", &controlPairToTs);
		tree->Branch("randomPairToTs", "std::vector<std::pair<double,double>>", &randomPairToTs);
		tree->Branch("randomIsolatedPairToTs", "std::vector<std::pair<double,double>>", &randomIsolatedPairToTs);
	}
	~SourceClusterWriter() {
		file->Write();
	}

	std::unique_ptr<TFile> file;
	TTree* tree;

	std::vector<PositionHit> hits{};
	int nHits{}, clusterChip{};
	double meanCol{}, meanRow{}, meanToT{}, meanToA{};
	double stdDevCol{}, stdDevRow{};
	Long64_t frameToA{};

	int nNeighbours{};
	std::vector<double> isolatedToTs;
	std::vector<std::pair<double,double> > isolatedPairToTs, controlPairToTs, randomPairToTs, randomIsolatedPairToTs;

};

std::vector<PositionHit> convertHitsQuad( std::vector<Hit>* chips [], double driftSpeed, double t0Offset=0 /*value to subtract from time*/) {
	std::vector<PositionHit> hits;
	for(int i=0; i<4; i++) {
		auto chipHits=convertHitsTPC(*chips[i], i, driftSpeed, t0Offset); //t0 is subtracted from hit time
		hits.insert(hits.end(), chipHits.begin(), chipHits.end() );
	}
	return hits;
}


template<class T > //=std::list<HitCluster>
inline TProfile2D profile2DPixel(const T& cluster) {
	TProfile2D prof("pixelProfile", "pixelProfile", 512,0,512,512,0,512);
	for(auto& iHit : cluster ) {
		auto h=iHit;
		prof.Fill(
				(h.chip>=2)*256+(h.chip==0||h.chip==3)*(256-2*h.column)+h.column,
				(h.chip==1||h.chip==2)+255+(h.chip==0||h.chip==3)*-2*h.row+h.row,
				h.ToT*0.025);
	}
	return prof;
}

void drawHitMapProfile(TProfile2D* axisObject) {
	gStyle->SetOptTitle(0);
	gStyle->SetPalette(kInvertedDarkBodyRadiator);
	double totAxis=1.6;
	static TCanvas* canv=new TCanvas("cluster_display2DPixel", "Event display 2D per pixel", 900,800);
	canv->cd();

	auto yaxis=axisObject->GetYaxis();
	if(!yaxis) {std::cout<<"could not get yaxis!?\n"; throw 1;}
	yaxis->SetTitle("Rows");
	auto xaxis=axisObject->GetXaxis();
	if(!xaxis) {std::cout<<"could not get xaxis!?\n"; throw 1;}
	xaxis->SetTitle("Columns");
	for(auto axis : {xaxis, yaxis} )
		axis->SetTitleOffset(1.1);
	for(auto axis : {xaxis, yaxis} ) {
		axis->SetTitleSize(0.05);
		axis->SetLabelSize(0.05);
	}
	axisObject->SetMaximum(totAxis);
	axisObject->SetMinimum(0);

	axisObject->DrawCopy("colz0");
	//pointTree.Draw( "h.position.y:h.position.x" , "", "same");
	gPad->SetMargin(0.15,0.1,0.1,0.05);//l r b t
	gPad->Update();

	gPad->Update();

}

template<typename T>
double getMean( std::vector<PositionHit> posHits, T Hit::*member) {
	return std::accumulate( posHits.begin(), posHits.end(), 0 , [&member](double x, const PositionHit& h){ return x+h.*member; } )/posHits.size();
}

std::vector<std::pair<double, double> > makeRandomPairs(
		std::vector<double>& ToTs) {
	//get random ToT pairs
	std::random_shuffle(ToTs.begin(), ToTs.end());
	std::vector<std::pair<double, double> > randomPairToTs;
	for (unsigned i = 0; i + 1 < ToTs.size(); i++) {
		randomPairToTs.push_back(std::minmax(ToTs[i], ToTs[i + 1]));
	}
	return randomPairToTs;
}



bool dropHit(std::string outputFileName) {
	const std::map<std::string, double> expectedNHits {
			{"run774_clusters_droph.root",224},
			{"run775_clusters_droph.root",275},
			{"run776_clusters_droph.root",359},
			{"run777_clusters_droph.root",506}
	};
	double target=220;
	if(expectedNHits.find(outputFileName)==expectedNHits.end()) return false;
	double expected=expectedNHits.at(outputFileName);
	if(expected<target) return false;
	return target/expected < gRandom->Uniform();
}

void processSource(std::string fileName, std::string outputFileName="output.root") {

	QuadTreeReader reader("data", fileName);
	SourceClusterWriter writer(outputFileName);

	Long64_t triggerT0=0;

	std::vector<PositionHit> posHits;
	while(reader.getNext() and reader.entryNumber<1E4) {
		auto newHits=convertHitsQuad(reader.chip, 0.055);

		auto tDiff=(reader.triggerToA-writer.frameToA)/4096.*25E-3;
		auto meanRow=getMean(newHits, &Hit::row);
		auto meanCol=getMean(newHits, &Hit::column);

		if( //another frame of same signal
					fabs( tDiff - 409.6)<10 //frame difference of 409.6 us
					and fabs( meanRow - writer.meanRow)<10
					and fabs( meanCol - writer.meanCol)<10) {
//				std::cout<<"append to frame!\n";
				for(auto& h:newHits) h.nShiftedTrigger=posHits.back().nShiftedTrigger+1;
				posHits.insert(posHits.end(), newHits.begin(), newHits.end());
		} else {//new entry
//			std::cout<<"start new frame!\n";
			writer.tree->Fill();
			posHits=newHits;
		}

		for(auto& h : posHits) if(h.ToT*.025<0.1) h.flag=PositionHit::Flag::lowToT;
		auto nHitsPerChip=countHitsPerChip(posHits, true);

		//reject event if every chip has many hits
		if( std::count_if(nHitsPerChip.begin(), nHitsPerChip.end(), [](int i){return i>10;} )>1 ) continue;
		if( *std::max_element(nHitsPerChip.begin(), nHitsPerChip.end()) <30) continue;
		int chipWithCluster=std::max_element(nHitsPerChip.begin(), nHitsPerChip.end())-nHitsPerChip.begin();

		auto prof=profile2DPixel(posHits);
		if(prof.GetStdDev(1)>30 || prof.GetStdDev(2)>30) continue;

//		drawHitMapProfile(&prof);
//		if(not triggerT0) triggerT0=reader.triggerToA;
//		if(std::cin.get()=='q') break;




		writer.hits=posHits;
		writer.nHits=posHits.size();
		writer.clusterChip=chipWithCluster;
		writer.meanCol=getMean(newHits, &Hit::column);
		writer.meanRow=getMean(newHits, &Hit::row);
		writer.meanToT=getMean(newHits, &Hit::ToT);
		writer.meanToA=std::accumulate( posHits.begin(), posHits.end(), 0 , [](double x, const PositionHit& h){ return x+h.driftTime; } )/posHits.size();
		writer.stdDevCol=prof.GetStdDev(1);
		writer.stdDevRow=prof.GetStdDev(2);
		writer.frameToA=reader.triggerToA;

		const bool calculateCrossTalk=true;
		if(calculateCrossTalk) {
			CrossTalkCalculator crossTalk;
			std::vector<double> ToTs;
			int nHits=0;
			for(const auto& h : posHits) if(h.chip==chipWithCluster) {
				if( dropHit(outputFileName) ) continue;
				nHits++;
				crossTalk.fill(h.column, h.row, h.ToT);
				ToTs.push_back(h.ToT);
			}
			auto nNeighbours=crossTalk.countNeighbours();
			auto isolatedToTs=crossTalk.getIsolatedToTs();

			writer.nHits=nHits;//number of hits after dropping
			writer.nNeighbours=nNeighbours;
			writer.isolatedToTs=crossTalk.getIsolatedToTs();
			writer.isolatedPairToTs=crossTalk.getPairToTs(crossTalk.pair);
			writer.controlPairToTs=crossTalk.getPairToTs(crossTalk.controlPair);
			writer.randomPairToTs=makeRandomPairs(ToTs);
			writer.randomIsolatedPairToTs=makeRandomPairs(isolatedToTs);
		}
	}
}
