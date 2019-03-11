/*
 * combineTracks.h
 *
 *  Created on: Oct 15, 2018
 *      Author: cligtenb
 */

#ifndef TRACKCOMBINER_COMBINETRACKS_H_
#define TRACKCOMBINER_COMBINETRACKS_H_

#include "TrackCombiner.cpp"

void combineTracks(std::string quad, std::string mimosa, bool draw=false) {
	Alignment align("align.dat");
	TrackCombiner tc(quad, mimosa, align);
	tc.drawEvent=draw;

	tc.telescopeFitter.getEntry(10);
	auto hits=tc.telescopeFitter.getSpaceHits();

	if(hits.size() ) for(auto& h : hits[0]) std::cout<<h.ToT;

	tc.openFile(draw ? "tmp.root" : "combinedFit.root");
	tc.Process();
}

void updateAlignmentFromFile(std::string resultFile="combinedFit.root", std::string alignFile="align.dat") {
	Alignment alignment(alignFile);
	TFile file(resultFile.c_str(), "READ");
//		alignment.updateAll(file);
//	alignment.updateShifts(file);
		alignment.quad.updateShift(file,"quad/global");
//	alignment.updateRotations(file);
		alignment.quad.updateRotation(file,"quad");
//		alignment.updateTimeWalk(file);
//		alignment.updateDriftSpeed(file);
		alignment.write("align.dat");
}

void combineTracksAndUpdateAlignment(std::string quad, std::string mimosa, int nIterations=3) {
	for(int i=0; i<nIterations; i++) {
		combineTracks(quad, mimosa);
		updateAlignmentFromFile();
	}
}



#endif /* TRACKCOMBINER_COMBINETRACKS_H_ */
