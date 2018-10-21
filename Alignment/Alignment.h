/*
 * Alignment.h
 *
 *  Created on: Jul 18, 2018
 *      Author: cligtenb
 */

#ifndef LASERDATAFITTER_ALIGNMENT_H_
#define LASERDATAFITTER_ALIGNMENT_H_

#include "TGraph.h"
#include "TH1.h"
#include "TFitResult.h"


#include "AlignmentHolder.h"
#include "TimeWalkCorrector.h"
//#include "/user/cligtenb/rootmacros/getHistFromTree.h"
#include "../../rootmacros/getHistFromTree.h"


struct Alignment {
	std::array<ChipAlignment,4> chips{{0,1,2,3}};
	ShiftAndRotateAlignment quad{"QUAD"};
	TimeWalkCorrector timeWalk;
	AlignValueHolder<double> driftSpeed{"DRIFTSPEED"};
	HitErrorCalculator hitErrors;
	TelescopeAlignment mimosa;

	Alignment(std::string filename) {
		read(filename);
	}
	void read(std::string filename) {
		std::ifstream inFile(filename);
		if(not inFile) std::cerr<<"Could not open "<<filename<<"\n";
		for(auto& c : chips) c.read(inFile);
		quad.read(inFile);
		timeWalk.read(inFile);
		driftSpeed.read(inFile);
		hitErrors.read(inFile);
		mimosa.read(inFile);
	}
	void write(std::string filename) {
		std::ofstream outFile(filename);
		if(not outFile) std::cerr<<"Could not open "<<filename<<"\n";
		outFile.precision(10);
		for(auto& c : chips) c.write(outFile);
		quad.write(outFile);
		timeWalk.write(outFile);
		driftSpeed.write(outFile);
		hitErrors.write(outFile);
		mimosa.write(outFile);
	}
	void updateAll(TFile& ) ;
	void updateShifts(TFile& file);
	void updateRotations(TFile& file);
	void updateDriftSpeed(TFile& file);

	int getChipNumber(const TVector3& pos) const;
	TVector3 transformBack(const TVector3& pos) const;
	PositionHit& transform(PositionHit& pos) const;
	std::array<TVector3,4> getChipCorners(int chipIndex) const {
		auto corners=chips[chipIndex].getChipCorners();
		for(auto& c : corners) c=quad.rotateAndShift(c);
		return corners;
	}
	std::vector<TVector3> getAllChipCorners() const {
		std::vector<TVector3> all;
		for(int i=0; i<4; i++) {
			auto corners=getChipCorners(i);
			all.insert(all.end(), corners.begin(), corners.end());
		}
		return all;
	}
};

void Alignment::updateShifts(TFile& file) {
	quad.updateShift(file, "quad");
//	for (int i = 0; i < 4; i++) {
//		chips[i].updateShift(file, "chip" + std::to_string(i + 1));
//	}
}
void Alignment::updateRotations(TFile& file) {
	quad.updateRotation(file, "quad");
//	for (int i = 0; i < 4; i++) {
//		chips[i].updateRotation(file, "chip" + std::to_string(i + 1));
//	}
}

void Alignment::updateDriftSpeed(TFile& file) {
	//calculate driftspeed
	TTree* tree = getObjectFromFile<TTree>("fitResults", &file);
	TH1* hist = getHistFromTree(tree, "hitAverage.z:laser.z",
			"fabs(hitAverage.x)+fabs(hitAverage.y)>0", "driftprof(26,-0.5,25.5)", "prof goff");
	if (hist->GetStdDev() > 5) {
		auto fitResult = hist->Fit("pol1", "QS");
		std::cout << "update drift speed by a factor " << fitResult->Parameter(1)
				<< "\n";
		driftSpeed.value /= fitResult->Parameter(1);
	} else {
		std::cout << "skipped updating drift speed\n";
	}
}

void Alignment::updateAll(TFile& file) {
	updateShifts(file);
	updateRotations(file);

	timeWalk.update(file);

	updateDriftSpeed(file);
}

inline int Alignment::getChipNumber(const TVector3& pos) const {
	for(unsigned i=0; i<chips.size(); i++) {
		if(chips[i].isOnChip(pos)) return i+1;
	}
	return 0;
}

inline TVector3 Alignment::transformBack(const TVector3& inPos) const {
	auto pos=quad.rotateAndShiftBack(inPos);
	auto chipNumber=getChipNumber(pos);
//	std::cout<<" chip "<<chipNumber<<": ";
	if(chipNumber) {
		pos=chips[chipNumber-1].rotateAndShiftBack(pos);
	} else {
		pos=TVector3{0,0,0};
	}
	return pos;
}


inline PositionHit& Alignment::transform(PositionHit& h) const {
	h.position=chips[h.chip].rotateAndShift(h.position);
	h.position=quad.rotateAndShift(h.position);
	return h;
}

#endif /* LASERDATAFITTER_ALIGNMENT_H_ */
