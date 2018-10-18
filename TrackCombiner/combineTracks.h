/*
 * combineTracks.h
 *
 *  Created on: Oct 15, 2018
 *      Author: cligtenb
 */

#ifndef TRACKCOMBINER_COMBINETRACKS_H_
#define TRACKCOMBINER_COMBINETRACKS_H_

#include "TrackCombiner.cpp"

void combineTracks(std::string quad, std::string mimosa) {
	Alignment align("align.dat");
	TrackCombiner tc(quad, mimosa, align);

	tc.telescopeFitter.getEntry(10);
	auto hits=tc.telescopeFitter.getSpaceHits();

	if(hits.size() ) for(auto& h : hits[0]) std::cout<<h.ToT;

	tc.openFile("combinedFit.root");
	tc.Process();
}



#endif /* TRACKCOMBINER_COMBINETRACKS_H_ */
