/*
 * combineTracks.h
 *
 *  Created on: Oct 15, 2018
 *      Author: cligtenb
 */

#ifndef TRACKCOMBINER_COMBINETRACKS_H_
#define TRACKCOMBINER_COMBINETRACKS_H_

#include "TrackCombiner.cpp"

void combineTracks() {
	Alignment align("align.dat");
	TrackCombiner tc("../eventBuilder/run668_events.root", "20181005-203345_M26_TELESCOPE_combined.root", align);
	tc.openFile("combinedFit.root");
	tc.Process();
}



#endif /* TRACKCOMBINER_COMBINETRACKS_H_ */
