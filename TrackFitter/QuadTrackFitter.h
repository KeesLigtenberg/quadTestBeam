/*
 * LaserDataFitter.h
 *
 *  Created on: Jul 16, 2018
 *      Author: cligtenb
 */

#ifndef LASERDATAFITTER_LASERDATAFITTER_H_
#define LASERDATAFITTER_LASERDATAFITTER_H_

#include <string>
#include <vector>
#include <memory>

#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TVector3.h"

#include "../eventBuilder/Hit.h"
#include "PositionHit.h"
#include "../Alignment/Alignment.h"

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

	std::unique_ptr<TFile> file;
	TTree* tree;

	Long64_t triggerToA{};
  unsigned triggerNumber{};
	std::vector<Hit>* chip[4]{};
};

struct ChipHistogrammer;

class QuadTrackFitter {
public:
	QuadTrackFitter(std::string fileName);
	virtual ~QuadTrackFitter();
	void Loop(std::string outputFile, const Alignment& alignment);
	int getEntry(Long64_t entryNumber) {
		if(entryNumber>=reader.tree->GetEntries()) {
			return false;
			std::cout<<"\nReached end of quad tree at entry "<<entryNumber<<"\n";
		}
		auto retStatus= reader.tree->GetEntry(entryNumber);
		if(!retStatus) std::cout<<"\ncould not get entry "<<entryNumber<<" from quad tree\n";

		return retStatus;
	} //returns true if valid

	QuadTreeReader reader;
	std::vector<PositionHit> posHits;

	unsigned triggerNumber() { return (reader.triggerNumber); }
	Long64_t numberOfEntries() { return reader.tree->GetEntries(); }


	std::vector<PositionHit> getSpaceHits(const Alignment& alignment);
};




#endif /* LASERDATAFITTER_LASERDATAFITTER_H_ */
