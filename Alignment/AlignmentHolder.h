/*
 * AlignmentHolder.h
 *
 *  Created on: Jul 23, 2018
 *      Author: cligtenb
 */

#ifndef LASERDATAFITTER_ALIGNMENTHOLDER_H_
#define LASERDATAFITTER_ALIGNMENTHOLDER_H_


#include <string>
#include <array>
#include <iostream>

#include "TVector3.h"
#include "TH1.h"
#include "TH2.h"

#include "../TrackFitter/PositionHit.h"

//#include "/user/cligtenb/rootmacros/getObjectFromFile.h"//
#include "../../rootmacros/getObjectFromFile.h"
#include "../TelescopeTrackFitter/getMeanFromGausFit.h"


//should actually be in namespace
	std::istream& operator>>(std::istream& in, TVector3& v) {
		double x,y,z;
		in>>x>>y>>z;
		v.SetXYZ(x,y,z);
		return in;
	}
	std::ostream& operator<<(std::ostream& out, const TVector3& v) {
		out<<v.x()<<" "<<v.y()<<" "<<v.z();
		return out;
	}

namespace { //todo:move to cpp file
	//line segment from line1 to line2, distance to p
	double getDistanceToLineSegment(TVector2 line1, TVector2 line2, TVector2 p) {
		line2-=line1;
		p-=line1;
		double projected=line2.Unit()*p;
		if(projected<0) return p.Mod();
		if(projected>line2.Mod() ) return (p-line2).Mod();
		return (p-projected*line2.Unit()).Mod();
	}

	//intersect line segments l1, l2, m1, m2
	TVector3 orientation(TVector3 q, TVector3 r, TVector3 p) {
		return (q-r).Cross(p-r);
	}
	//segments are expected to lie in the same plane
	bool doesIntersect( TVector3 l1, TVector3 l2, TVector3 m1, TVector3 m2 ){
		bool intersects=
				orientation(l2,l1,m1)*orientation(l2,l1,m2)<0 //orientation of (line x m1) and (line x m2) is different
				and
				orientation(m2,m1,l1)*orientation(m2,m1,l2)<0; //orientation of (mline x l1) and (mline x l2) is different
//		std::cout<<"doesIntersect("<<l1<<"; "<<l2<<"; "<<m1<<"; "<<m2<<")="<<intersects<<"\n";
		return intersects;
	}

}

struct AlignmentHolder {
	AlignmentHolder(std::string name) : name(name) {}
	virtual ~AlignmentHolder() {}

	void read(std::istream&);
	void write(std::ostream&);

private:
	virtual void readParameters(std::istream&) =0;
	virtual void writeParameters(std::ostream&) =0;

	const std::string name;
};

void AlignmentHolder::read(std::istream& in) {
	//check header
	std::string word;
	if(not in or not (in>>word) ) std::cerr<<"could not read word\n";
	if(word!=name) std::cerr<<"wrong header: expected "<<name<<", but found "<<word<<"\n";
	//continue to read the parameters
	readParameters(in);
}
void AlignmentHolder::write(std::ostream& out) {
	out<<name<<"\n";
	writeParameters(out);
}

template<class T>
struct AlignValueHolder : AlignmentHolder {
	using AlignmentHolder::AlignmentHolder;//use Alignmentholder constructor
	T value;
	void readParameters(std::istream& in) {
		in>>value;
	}
	void writeParameters(std::ostream& out) {
		out<<value<<"\n";
	}
};

struct ShiftAndRotateAlignment : AlignmentHolder {
	using AlignmentHolder::AlignmentHolder;
	TVector3 shift{};
	TVector3 COM{};
	TVector3 rotation{};

	virtual Vec3 rotateAndShift(const Vec3& hpos) const {
		TVector3 Thpos=hpos; //convert to TVector3 and then call the same function
		rotateAndShift(Thpos);
		return Thpos;
	}
	virtual TVector3& rotateAndShift(TVector3& hpos) const {
		hpos+=shift;
		for(int i=0; i<3; i++) {
			const std::array<TVector3, 3> unitVectors={ TVector3{1,0,0}, TVector3{0,1,0}, TVector3{0,0,1} };
			hpos=RotateAroundPoint(hpos, rotation[i],COM, unitVectors[i]);
		}
		return hpos;
	}
	virtual Vec3 rotateAndShiftBack(const Vec3& hpos) const {
		TVector3 Thpos=hpos; //convert to TVector3 and then call the same function
		rotateAndShiftBack(Thpos);
		return Thpos;
	}
	virtual TVector3& rotateAndShiftBack(TVector3& hpos) const {
		for(int i=2; i>=0; i--) {
			const std::array<TVector3, 3> unitVectors={ TVector3{1,0,0}, TVector3{0,1,0}, TVector3{0,0,1} };
			hpos=RotateAroundPoint(hpos, -rotation[i],COM, unitVectors[i]);
		}
		hpos-=shift;
		return hpos;
	}

	void updateShift(TFile&, std::string dirName);
	void updateShift(TFile& file, const std::string& dirName, int axis);
	void updateRotation(TFile&, std::string dirName);

private:
	void readParameters(std::istream& in) {
		in>>shift>>COM>>rotation;
	}
	void writeParameters(std::ostream& out) {
		out<<shift<<"\n"<<COM<<"\n"<<rotation<<"\n";
	}
};

struct ChipAlignment : ShiftAndRotateAlignment {
	ChipAlignment(int chipNumber) : ShiftAndRotateAlignment("CHIP" + std::to_string(chipNumber) ), chipNumber(chipNumber) {}
	int chipNumber;//starting at 0

	using ShiftAndRotateAlignment::rotateAndShift;
	TVector3& rotateAndShift(TVector3& hpos) const {
		auto dir=(chipNumber==0 || chipNumber==3) ? -1 : 1 ; //because of chip orientation
		hpos[0]*=dir; hpos[1]*=dir;
		hpos+=shift;
		for(int i=0; i<3; i++) {
			const std::array<TVector3, 3> unitVectors={ TVector3{1,0,0}, TVector3{0,1,0}, TVector3{0,0,1} };
			hpos=RotateAroundPoint(hpos, rotation[i],COM, unitVectors[i]);
		}
		return hpos;
	}
	using ShiftAndRotateAlignment::rotateAndShiftBack;
	TVector3& rotateAndShiftBack(TVector3& hpos) const {
		for(int i=2; i>=0; i--) {
			const std::array<TVector3, 3> unitVectors={ TVector3{1,0,0}, TVector3{0,1,0}, TVector3{0,0,1} };
			hpos=RotateAroundPoint(hpos, -rotation[i],COM, unitVectors[i]);
		}
		hpos-=shift;
		auto dir=(chipNumber==0 || chipNumber==3) ? -1 : 1 ;//because of chip orientation
		hpos[0]*=-dir; hpos[1]*=dir;
		return hpos;
	}

	std::array<TVector3,4> getChipCorners() const { //fixme
		std::array<TVector3,4> corners={//upper-left, upper-right, bottom-left, bottom-right
				TVector3{0,0,0},
				TVector3{14.08,0,0},
				TVector3{14.08,14.08,0},
				TVector3{0,14.08,0}
		};
		for(auto& c: corners) {
			rotateAndShift(c);
		}
		return corners;
	}
	double getDistanceFromEdge(TVector3 pos) const {
				auto corners=getChipCorners();
				double distance=getDistanceToLineSegment(corners[0].XYvector(), corners[3].XYvector(), pos.XYvector());
				for(unsigned i=0; i<corners.size()-1; i++) {
					distance=std::min(distance, getDistanceToLineSegment(corners[i].XYvector(), corners[i+1].XYvector(), pos.XYvector()));
				}
//				std::cout<<distance<<"\n";
				return distance;
	}
	bool isOnChip(TVector3 pos) const {
		//if pos has 1 crossing with line segment to infinite it is inside, if 0 or 2 then not inside
		auto corners=getChipCorners();
		for(auto& c:corners) c.SetZ(0);
		TVector3 inf{0,0,0};
		pos.SetZ(0);
		int nCrossings=0;
		nCrossings+=doesIntersect(corners[0], corners[3], pos,inf);
		for(unsigned i=0; i<corners.size()-1; i++) {
			nCrossings+=doesIntersect(corners[i], corners[i+1], pos,inf);
		}
		return nCrossings==1;
	}

};

void ShiftAndRotateAlignment::updateShift(TFile& file,const std::string& dirName, int i) {
		const std::string histNames[] = { "xResidual", "yResidual", "zResidual" };
		auto hist = getObjectFromFile<TH1>(dirName + "/" + histNames[i], &file);
		const int minEntries = 100;
		if (hist->GetEntries() < minEntries) {
			std::cout << "skipped " << dirName << " because less than " << minEntries
					<< " in histogram\n";
		} else {
			double mean = hist->GetMean(); //getMeanFromGausCoreFit(*hist);
			shift[i] -= mean;
			std::cout << "update " << histNames[i] << " by " << mean << "\n";
		}
}

void ShiftAndRotateAlignment::updateShift(TFile& file, std::string dirName) {
	std::cout<<dirName<<"\n";
	for (int i = 0; i < 3; i++) {
		updateShift(file, dirName, i);
	}
}

void ShiftAndRotateAlignment::updateRotation(TFile& file, std::string dirName) {
	//Find Center of Mass (point of rotation)
	auto hitmap=getObjectFromFile<TH2>(dirName+"/positionHitMap", &file);
	auto zHit=getObjectFromFile<TH1>(dirName+"/zHit", &file);
	const int minEntries=10000;
	if(hitmap->GetEntries()<minEntries) {std::cout<<"skipped "<<dirName<<" because less than "<<minEntries<<" in histogram\n"; return;}
	double meanX=hitmap->GetMean(1), meanY=hitmap->GetMean(2), meanZ=zHit->GetMean();
	auto oldCOM=COM;
	COM.SetXYZ(meanX,meanY,meanZ);
	std::cout<<"updated COM by "<<COM-oldCOM<<"\n";

	//Find rotation
	for(int i=0; i<3; i++) {
		const std::array<std::string,3> x={"x", "y", "z"};
		auto xHist=getObjectFromFile<TH1>(dirName+"/"+x[i]+"Rotation", &file);
		std::cout<<"update "<<x[i]<<" rotation by "<<xHist->GetMean()<<"\n";
		rotation[i]+=xHist->GetMean();
	}

}

struct HitErrorCalculator : AlignmentHolder {
	HitErrorCalculator() : AlignmentHolder("HITERROR") {}
	TVector3 sigma0{1,1,1}, diffusion{};
	double z0{0};

	TVector3 hitError(double z) const {
		TVector3 error{};
		if(z<z0) z=z0;
		for(int i=0; i<3; i++) {
			error[i]=sqrt(sigma0[i]*sigma0[i]+diffusion[i]*(z-z0));
		}
		return error;
	}

private:
	void readParameters(std::istream& in) {
		in>>z0>>sigma0>>diffusion;
	}
	void writeParameters(std::ostream& out) {
		out<<z0<<"\n"<<sigma0<<"\n"<<diffusion<<"\n";
	}

};

struct TelescopeAlignment : AlignmentHolder {
	TelescopeAlignment() : AlignmentHolder("TELESCOPE") {};
	std::vector<std::pair<double,double>> shifts;
	std::vector<std::pair<double,double>> COMs;
	std::vector<double> angles;
	std::pair<double,double> slopes;

	const int nPlanes=6;

	virtual void readParameters(std::istream& in) {
		for(int i=0; i<nPlanes; i++) {
			std::string buf{};
			in>>buf;
			if(buf!="PLANE"+std::to_string(i)) throw 1;
			double shiftx{}, shifty{}, comx{}, comy{}, anglez{};
			in>>shiftx>>shifty>>comx>>comy>>anglez;
			shifts.push_back({shiftx,shifty});
			COMs.push_back({comx,comy});
			angles.push_back(anglez);
		}
		std::string buf{};
		in>>buf;
		if(buf!="SLOPES") throw 1;
		in>>slopes.first>>slopes.second;
	}
	virtual void writeParameters(std::ostream& out) {
		for(int i=0; i<nPlanes; i++) {
			out<<"PLANE"<<i<<"\n";
			out<<shifts[i].first<<" "<<shifts[i].second<<"\n";
			out<<COMs[i].first<<" "<<COMs[i].second<<"\n";
			out<<angles[i]<<"\n";
		}
		out<<"SLOPES\n";
		out<<slopes.first<<" "<<slopes.second<<"\n";
	}


};


#endif /* LASERDATAFITTER_ALIGNMENTHOLDER_H_ */
