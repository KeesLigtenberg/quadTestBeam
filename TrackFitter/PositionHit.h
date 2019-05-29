#ifndef POSITIONHIT_H
#define POSITIONHIT_H
//based on single chip TrackFitter/PositionHit but changes were made!

#include "../eventBuilder/Hit.h"

#include <algorithm>
#include <functional>
#include <set>
#include <memory>

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
	Vec3 residual{}, error{1,1,1};

	enum class Flag : int {valid=1, highResidualxy=-1, highResidualz=-3, lowToT=-2, smallz=-4, debug=-5, shiftedTrigger=-6, outsideFiducial=-7} flag=Flag::valid;
	bool isValid() const { return flag==Flag::valid; }


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
Container convertHitsTPC(const std::vector<Hit>& hv, unsigned char chip, double driftSpeed, double t0Offset, double pixelwidth=0.055, double pixelheight=0.055 ) {
		Container phv;
		for(auto h : hv) {
			// .5 to center position in pixel, removed?
			h.driftTime-=t0Offset;
			phv.emplace_back(
					PositionHit{ (h.column)*pixelwidth/*x*/, (h.row)*pixelheight /*y*/,h.driftTime*driftSpeed /*z*/, chip, h }
			);
		}
		return phv;
}


template< class Container=std::vector<PositionHit> >
Container convertHits(const std::vector<Hit>& hv, unsigned char chip, double planePosition, double pixelwidth=0.055, double pixelheight=0.055 ) {
		Container phv;
		for(auto& h : hv) {
			// .5 to center position in pixel, removed?
			phv.emplace_back(
					PositionHit{ (h.column)*pixelwidth/*x*/, (h.row)*pixelheight /*y*/,planePosition /*z*/, chip, h }
			);
		}
		return phv;
}

PositionHit& flagResidual(PositionHit& h, const TVector3& maxResidual) {
	if(fabs(h.residual.x)>maxResidual.x()
			or fabs(h.residual.y)>maxResidual.y() ) {
		h.flag=PositionHit::Flag::highResidualxy;
	}
	if(fabs(h.residual.z)>maxResidual.z()) {
		h.flag=PositionHit::Flag::highResidualz;
	}
	return h;
}

PositionHit& flagResidualPull(PositionHit& h, const TVector3& maxResidualPull) {
	for(int i=0; i<3; i++) {
		if(fabs(h.residual[i])/h.error[i]>maxResidualPull[i]) {
			h.flag= i==2 ? PositionHit::Flag::highResidualz : PositionHit::Flag::highResidualxy;
		}
	}
	return h;
}

PositionHit& flagShiftedTrigger(PositionHit& h, int maxShift=0) {
	if(h.nShiftedTrigger>maxShift) {
		h.flag = PositionHit::Flag::shiftedTrigger;
	}
	return h;
}

template<class Container>
TVector3 getAveragePositionCondition(Container hv, std::function<bool(const PositionHit&)> isValid) {
	TVector3 sum(0,0,0);
	int n=0;
	for(const PositionHit& hit : hv) {
		if(not isValid(hit)) continue;
		sum+=hit.position;
		n++;
	}
	sum*=(1./n);
	return sum;
}
template<class Container>
TVector3 getWeightedAveragePosition(Container hv, std::function<bool(const PositionHit&)> isValid) {
	TVector3 sum(0,0,0);
	TVector3 weightsum(0,0,0);
	for(const PositionHit& hit : hv) {
		if(not isValid(hit)) continue;
		for(int i=0; i<3; i++) {
			sum[i]+=hit.position[i]/hit.error[i]/hit.error[i];
			weightsum[i]+=1./hit.error[i]/hit.error[i];
		}
	}
	for(int i=0; i<3; i++)
		sum[i]/=weightsum[i];
	return sum;
}


template<class Container>
TVector3 getAveragePosition(Container hv, bool rejectFlagged=false) {
	if(rejectFlagged) {
		return getAveragePositionCondition(hv, [](const PositionHit& h){return h.flag==PositionHit::Flag::valid;} );
	} else {
		return getAveragePositionCondition(hv, [](const PositionHit& h){return true;});
	}
}

template<class Container>
TVector3 getAveragePosition(Container hv, std::set<PositionHit::Flag> rejectFlags) {
	return getAveragePositionCondition(hv, [&rejectFlags](const PositionHit& h){return rejectFlags.find(h.flag)==rejectFlags.end();});
}

template<class Container>
std::vector<Vec3> getWeightedAveragePositionPerChip(Container hv) {
	std::vector<Vec3> averages;
	for(int i=0; i<4; i++) {
		averages.push_back( getWeightedAveragePosition(hv, [i](const PositionHit& h){ return h.isValid() && h.chip == i; } ) );
	}
	return averages;
}

std::vector<int> countHitsPerChip(const std::vector<PositionHit>& hits, bool rejectFlagged=false) {
	std::vector<int> nHits(4);
	for(auto& h : hits) {
		if(rejectFlagged and h.flag!=PositionHit::Flag::valid) continue;
		nHits[h.chip]++;
	}
	return nHits;
}

std::vector<std::vector<PositionHit> > getHitsPerChip(const std::vector<PositionHit>& hits, bool rejectFlagged=true) {
	std::vector<std::vector<PositionHit>> perChip(4);
	for(const auto& h : hits) {
		if(rejectFlagged and not h.isValid()) continue;
		perChip[h.chip].push_back(h);
	}
	return perChip;
}

int countTotalValidHits(const std::vector<PositionHit>& hits) { //simply use size() to get to total number of unflagged hits!
	return std::count_if(hits.begin(), hits.end(), [](const PositionHit& h){return h.flag==PositionHit::Flag::valid; });
}

#pragma link C++ class PositionHit+;
#pragma link C++ class std::vector<PositionHit>+;
#pragma link C++ class Residual+;

#endif
