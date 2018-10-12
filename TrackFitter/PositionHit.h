#ifndef POSITIONHIT_H
#define POSITIONHIT_H
//based on single chip TrackFitter/PositionHit but changes were made!

#include "../eventBuilder/Hit.h"
#include "TVector3.h"

TVector3& RotateAroundPoint(TVector3& v, double rotation, const TVector3& rotationPoint, const TVector3& rotationAxis) {
	v-=rotationPoint;
	v.Rotate(rotation, rotationAxis);
	v+=rotationPoint;
	return v;
};

struct PositionHit : Hit{
	PositionHit() : PositionHit(0,0,0,0) {};
	PositionHit(double x, double y, double z, unsigned char chip, const Hit& hit=Hit()) :
			position(x,y,z), chip(chip), Hit(hit) {};
	Vec3 position;
	unsigned char chip;
	Vec3 residual{}, error{};

	enum class Flag : int {valid=1, highResidualxy=-1, highResidualz=-3, lowToT=-2, smallz=-4} flag=Flag::valid;


	void RotatePosition(double rotation, const TVector3& rotationPoint, const TVector3& rotationAxis) {
		TVector3 v=position;
		v=RotateAroundPoint(v, rotation, rotationPoint, rotationAxis);
		position=v;
	};
	const Vec3& calculateResidual(const TVector3& expectedPosition) {
		residual=position-expectedPosition;
		return residual;
	}
};

template< class Container=std::vector<PositionHit> >
Container convertHits(const std::vector<Hit>& hv, unsigned char chip, double driftSpeed, double pixelwidth=0.055, double pixelheight=0.055 ) {
		Container phv;
		for(auto& h : hv) {
			// .5 to center position in pixel, removed?
			phv.emplace_back(
					PositionHit{ (h.column)*pixelwidth/*x*/, (h.row)*pixelheight /*y*/,h.driftTime*driftSpeed /*z*/, chip, h }
			);
		}
		return phv;
}


#pragma link C++ class PositionHit+;
#pragma link C++ class std::vector<PositionHit>+;
#pragma link C++ class Residual+;

#endif
