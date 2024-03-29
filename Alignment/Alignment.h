/*
 * Alignment.h
 *
 *  Created on: Jul 18, 2018
 *      Author: cligtenb
 */

#ifndef LASERDATAFITTER_ALIGNMENT_H_
#define LASERDATAFITTER_ALIGNMENT_H_

#include <array>
#include <set>
#include <iostream>
#include <fstream>

#include "TGraph.h"
#include "TH1.h"
#include "TFitResult.h"
#include "TPolyLine.h"


#include "AlignmentHolder.h"
#include "TimeWalkCorrector.h"
#include "/user/cligtenb/rootmacros/getHistFromTree.h"
//#include "../../rootmacros/getHistFromTree.h"


struct Alignment {
	std::array<ChipAlignment,4> chips{{0,1,2,3}};
	ShiftAndRotateAlignment quad{"QUAD"};
	std::array<TimeWalkCorrector,4> timeWalks{ {"TIMEWALKPARAMETERS0", "TIMEWALKPARAMETERS1", "TIMEWALKPARAMETERS2", "TIMEWALKPARAMETERS3"} };
	AlignValueHolder<double> driftSpeed{"DRIFTSPEED"};
	AlignValueHolder<double> t0Offset{"T0OFFSET"};
	HitErrorCalculator hitErrors;
	TelescopeAlignment mimosa;
	std::array<ToTCorrection,4> totCorrections{ {"TOTCORRECTION0","TOTCORRECTION1","TOTCORRECTION2","TOTCORRECTION3"} };

	Alignment(std::string filename) {
		read(filename);
	}
	void read(std::string filename) {
		std::ifstream inFile(filename);
		if(not inFile) std::cerr<<"Could not open "<<filename<<"\n";
		for(auto& c : chips) c.read(inFile);
		quad.read(inFile);
		for(auto& t: timeWalks) t.read(inFile);
		driftSpeed.read(inFile);
		t0Offset.read(inFile);
		hitErrors.read(inFile);
		mimosa.read(inFile);
		for(auto& t : totCorrections) t.read(inFile);
	}
	void write(std::string filename) {
		std::ofstream outFile(filename);
		if(not outFile) std::cerr<<"Could not open "<<filename<<"\n";
		outFile.precision(10);
		for(auto& c : chips) c.write(outFile);
		quad.write(outFile);
		for(auto& t: timeWalks) t.write(outFile);
		driftSpeed.write(outFile);
		t0Offset.write(outFile);
		hitErrors.write(outFile);
		mimosa.write(outFile);
		for(auto& t : totCorrections) t.write(outFile);
	}
	void updateAll(TFile& ) ;
	void updateShifts(TFile& file);
	void updateRotations(TFile& file);
	void updateShears(TFile& file);
	void updateDriftSpeed(TFile& file);
	void updateTimeWalk(TFile& file);

	void resetQuadOrigin();

	int getChipNumber(const TVector3& pos) const;
	TVector3 transformBack(const TVector3& pos) const;
	PositionHit& transform(PositionHit& pos) const;
	std::array<TVector3,4> getChipCorners(int chipIndex) const {
		auto corners=chips[chipIndex].getChipCorners();
		for(auto& c : corners) c=quad.rotateAndShift(c);
		return corners;
	}
	std::vector<TVector3> getAllChipCorners() const { //in global frame
		std::vector<TVector3> all;
		for(int i=0; i<4; i++) {
			auto corners=getChipCorners(i);
			all.insert(all.end(), corners.begin(), corners.end());
		}
		return all;
	}
	std::vector<TVector3> getAllChipCornersQuad() const { //in quad frame
		std::vector<TVector3> all;
		for(int i=0; i<4; i++) {
			auto corners=chips[i].getChipCorners();
			all.insert(all.end(), corners.begin(), corners.end());
		}
		return all;
	}
	void drawChipEdges(bool globalFrame=true, Color_t color=kBlack) const { //in global frame
		for(int i=0; i<4; i++) {
			auto corners=globalFrame ? getChipCorners(i) : chips[i].getChipCorners();
			TPolyLine l;
			for(auto& corner : corners) {
				l.SetNextPoint(corner.x(), corner.y());
			}
			l.SetNextPoint(corners[0].x(), corners[0].y());
			l.SetLineColor(color);
			l.DrawClone();
		}
	}
};


inline double getDriftSpeedFactor(TFile& file, bool draw) {
	//calculate driftspeed
//	std::cout<<"get zResidualByz from quad\n";
	std::string histogramName="zResidualByz";
	TH2D* zResidualByz = getObjectFromFile<TH2D>("quad/locExp/"+histogramName, &file);

//	std::cout<<"fit slices!\n";
	TF1* gausRange=new TF1("gausRange","gaus(0)", -0.4,0.4);
	gausRange->SetParameters(4E4,0.05,0.22);
	zResidualByz->FitSlicesY(gausRange,0,-1, 5, "QNR");

//	std::cout<<"get slices!\n";
	TH1 *zResidualByz_1 = dynamic_cast<TH1*>( gDirectory->Get( (histogramName+"_1").c_str() ) );
	if(!zResidualByz_1) {
		std::cout<<"could not get fitted slices!";
		return 0.;
	}
//	zResidualByz_1->Draw(""); gPad->Update(); std::cin.get();

	std::cout<<"fit pol1 to slices!\n";
	auto projectionx=zResidualByz->ProjectionX();	//fit only to core, get from projection
	TF1* line=new TF1("line","pol1(0)", projectionx->GetMean()-projectionx->GetStdDev(), projectionx->GetMean()+projectionx->GetStdDev() );
	auto fitResult=zResidualByz_1->Fit(line, "SQ");
//	gPad->Update(); std::cin.get();
//	std::cout<<"fit result: "<<fitResult<<"\n";

	auto slope=fitResult->Parameter(1);
	return slope;
}


void Alignment::updateShifts(TFile& file) {
//	quad.updateShift(file, "quad");
	for (int i = 0; i < 4; i++) { //skip chip0, because otherwise we would have a redundant d.o.f.
		chips[i].updateShift(file, "chip" + std::to_string(i));
	}
}
void Alignment::updateRotations(TFile& file) {
//	quad.updateRotation(file, "quad");
	for (int i = 0; i < 4; i++) {
		chips[i].updateRotation(file, "chip" + std::to_string(i));
	}
}

void Alignment::updateShears(TFile& file) {
//	quad.updateRotation(file, "quad");
	for (int i = 0; i < 4; i++) {
		chips[i].updateShear(file, "chip" + std::to_string(i));
	}
}


void Alignment::updateDriftSpeed(TFile& file) {
	auto slope=getDriftSpeedFactor(file,false);
	driftSpeed.value *= (1-slope);
	std::cout << "update drift speed by a factor " << (1-slope) <<"\n";
}


void Alignment::updateTimeWalk(TFile& file) {
	for(unsigned i=0; i<timeWalks.size(); i++) {
		TH2D* uncorrected=getObjectFromFile<TH2D>("chip"+std::to_string(i)+"/zResidualByToT", &file);
		timeWalks[i].update(uncorrected);
	}
}

void Alignment::updateAll(TFile& file) {
	updateShifts(file);
	updateRotations(file);

	updateTimeWalk(file);

	updateDriftSpeed(file);
}

void Alignment::resetQuadOrigin() {
	//reset quad origin to the top left corner of chip0 (center of first pixel, last column)

	//in global frame, move quad origin to upperleft corner
	auto upperLeft = getChipCorners(0)[1];
	quad.shift=upperLeft;

	//in local frame, move all chips, such that chip0 top left is on origin
	auto chip0shift= chips[0].getChipCorners()[1];
	for(auto& c : chips) {
		c.shift-=chip0shift;
	}

}

//quad coordinates, returns 0 if no chip was found!
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
