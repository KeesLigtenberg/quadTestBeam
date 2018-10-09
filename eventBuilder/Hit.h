#ifndef HIT_H
#define HIT_H

#include <vector>
#include "TVector3.h"

//Note: slightly different from Hit in testBeamTrackFitter!
struct Hit {
	Hit() : row(0), column(0), ToT(0), driftTime(0) {};
	Hit(unsigned char chipNumber, unsigned short row, unsigned short column, unsigned short ToT, int driftTime) : chip(chipNumber), row(row), column(column), ToT(ToT), driftTime(driftTime) {};
	virtual ~Hit() {}
	unsigned short row, column, ToT; //charge=ToT
	int driftTime;
};

struct Vec3 {//because TVector3 has a bug in combination with TTreeReader
	Vec3() : x(0), y(0), z(0) {}
	Vec3(double x,double y,double z) : x(x), y(y), z(z) {}
	Vec3(const TVector3& vec) : x(vec.x()), y(vec.y()), z(vec.z()) {}
	operator TVector3() const {return TVector3(x,y,z);}
	double operator[] (int i) { return TVector3(*this)[i]; }
	double x,y,z;
};

#pragma link C++ class Hit+;
#pragma link C++ class Vec3+;
#pragma link C++ class std::vector<Hit>+;

#endif
