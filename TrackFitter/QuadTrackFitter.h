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
				reader(treeName.c_str(), file.get() ) {}
	std::unique_ptr<TFile> file;
	TTreeReader reader;
  TTreeReaderValue<Long64_t> triggerToA{reader, "triggerToA"};
	TTreeReaderValue<std::vector<Hit> > chip[4] = {
			{reader, "chip0"},
			{reader, "chip1"},
			{reader, "chip2"},
			{reader, "chip3"}
	};
};



class QuadTrackFitter {
public:
	QuadTrackFitter(std::string fileName);
	virtual ~QuadTrackFitter();
	void Loop(std::string outputFile, const Alignment& alignment);

	QuadTreeReader tree;
	std::vector<PositionHit> posHits;

};

#endif /* LASERDATAFITTER_LASERDATAFITTER_H_ */
