/*
 * fitLaserData.h
 *
 *  Created on: Jul 26, 2018
 *      Author: cligtenb
 */

#ifndef LASERDATAFITTER_FITLASERDATA_H_
#define LASERDATAFITTER_FITLASERDATA_H_

#include "QuadTrackFitter.cpp"

void fitTracks(std::string inputFileName="tree.root", std::string outputFileName="fitted.root", std::string alignFile="align.dat") {
	QuadTrackFitter ldf{inputFileName};
	Alignment alignment{alignFile};
	//do fit
	ldf.Loop(outputFileName, alignment);
}

void fitAndUpdateAlignment(std::string inputFileName="tree.root", std::string outputFileName="fitted.root", std::string alignFile="align.dat", int nRepeat=1) {

	for(int i=0; i<nRepeat; i++)  {
		std::cout<<"iteration "<<i<<"\n\n";

		QuadTrackFitter ldf{inputFileName};
		Alignment alignment{alignFile};

		//do fit
		ldf.Loop(outputFileName, alignment);

		TFile file(outputFileName.c_str(), "READ");
//		alignment.updateAll(file);
//		alignment.updateShifts(file);
//		alignment.quad.updateShift(file,"quad",2);
//		alignment.updateRotations(file);
//		alignment.timeWalk.update(file);
//		alignment.updateDriftSpeed(file);
//		alignment.write("align.dat");

	}

}


#endif /* LASERDATAFITTER_FITLASERDATA_H_ */
