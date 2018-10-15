/*
 * transformHits.h
 *
 *  Created on: Jun 16, 2017
 *      Author: cligtenb
 */

#ifndef TRANSFORMHITS_H_
#define TRANSFORMHITS_H_

#include "../TrackFitter/PositionHit.h"

inline PositionHit& translateHit(PositionHit& h, double dx, double dy ) {
	h.position.x-=dx;
	h.position.y-=dy;
	return h;
}
inline PositionHit& translateHit(PositionHit& h, const std::pair<double,double>& translation) { return translateHit(h,translation.first, translation.second); }

//rotate around line parellel to z-axis
inline PositionHit& rotateHit(PositionHit& h, double rotation, const std::pair<double, double>& rotationPoint) {
	double sinr=sin(rotation), cosr=cos(rotation);

	const double xc=rotationPoint.first, yc=rotationPoint.second; //x and y center

	h.position.x=cosr*(h.position.x-xc)-sinr*(h.position.y-yc)+xc;
	h.position.y=cosr*(h.position.y-yc)+sinr*(h.position.x-xc)+yc;
	return h;
}

template<class hitCollection>
hitCollection& translateHits( hitCollection& hits, const std::vector<std::pair<double,double> >& translation ) {
	for(PositionHit& h : hits) translateHit(h, translation[h.chip] );
	return hits;
};

template<class hitCollection>
hitCollection& translateHits( hitCollection& hits, const std::pair<double,double>& translation ) {
	for(PositionHit& h : hits) translateHit(h, translation );
	return hits;
};

template<>
std::vector<std::vector<PositionHit> >& translateHits( std::vector<std::vector<PositionHit> >& hits, const std::vector<std::pair<double,double> >& translation ) {
	for(auto& v: hits)
		for(PositionHit& h : v)
			h=translateHit(h, translation[h.chip] );
	return hits;
};

template<class hitCollection>
hitCollection& rotateHits(hitCollection& hits, const std::vector<double>& rotations, const std::vector<std::pair<double, double>>& rotationPoints ) {
	for(PositionHit& h : hits) {
		h=rotateHit(h,rotations[h.chip], rotationPoints[h.chip]);
	}
	return hits;
}

template<class hitCollection>
hitCollection& rotateHits(hitCollection& hits, double rotation, const std::pair<double, double>& rotationPoints ) {
	for(PositionHit& h : hits) {
		h=rotateHit(h,rotation, rotationPoints);
	}
	return hits;
}

template<>
std::vector<std::vector<PositionHit> >& rotateHits(std::vector<std::vector<PositionHit> >& hits, const std::vector<double>& rotations, const std::vector<std::pair<double, double>>& rotationPoints) {
	for(auto& v : hits) {
		v=rotateHits(v,rotations,rotationPoints);
	}
	return hits;
}

#endif /* TRANSFORMHITS_H_ */
