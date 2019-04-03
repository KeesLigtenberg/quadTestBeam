/*
 * linearRegressionFit.cpp
 *
 *  Created on: Jun 15, 2017
 *      Author: cligtenb
 */

#include <iostream>
#include "TPolyLine.h"

#include "linearRegressionFit.h"

ClassImp(FitResult3D);

void FitResult3D::draw(double zmin, double zmax, Color_t color) const { //in timepix frame, with telescope coordinates!
	const int npoints=2;
	double x[npoints] = { XZ.at(zmin), XZ.at(zmax)};
	double y[npoints] = { YZ.at(zmin), YZ.at(zmax)};
	double z[npoints] = { zmin, zmax};
//	std::cout<<"Line "<<color<<" :";
//	for(int i=0; i<npoints; i++) std::cout<<" ("<<x[i]<<", "<<y[i]<<", "<<z[i]<<")";
//	std::cout<<"\n";
	TPolyLine3D l( npoints, x, z, y );
	l.SetLineColor(color);
	l.SetLineWidth(2);
	l.DrawClone();
}
void FitResult2D::draw(double zmin, double zmax, Color_t color) const { //in timepix frame, with telescope coordinates!
	const int npoints=2;
	double x[npoints] = {at(zmin), at(zmax)};
	double z[npoints] = { zmin, zmax};
	TPolyLine l( npoints, x, z );
	l.SetLineColor(color);
	l.SetLineWidth(2);
	l.DrawClone();
}
FitResult3D FitResult3D::makeRotated(double rotation, const TVector3& rotationPoint, const TVector3& rotationAxis ) const {
	if( fabs(XZ.interceptz-YZ.interceptz) > 1E-10 ) {std::cerr<<"intercepts should be described at the same point in space!\n"; throw "fabs(XZ.interceptz-YZ.interceptz) > 1E-10";}
	//rotate
	TVector3 slope(XZ.slope, YZ.slope, 1), intercept(XZ.intercept, YZ.intercept, XZ.interceptz);
	slope.Rotate(rotation, rotationAxis);
	intercept=RotateAroundPoint(intercept, rotation, rotationPoint, rotationAxis);
	//return result
	return FitResult3D{
		FitResult2D{
			slope.x()/slope.z(),
			-intercept.z()*slope.x()/slope.z()+intercept.x(),
			XZ.error,//todo:propagate errors!
			0,//intercept.z()
		},
		FitResult2D{
			slope.y()/slope.z(),
			-intercept.z()*slope.y()/slope.z()+intercept.y(),
			YZ.error,//todo:propagate errors!
			0,//intercept.z()
		}
	};
}

template<class HitContainer=HoughTransformer::HitCluster>
FitResult2D regressionXZ(const HitContainer& cluster, double interceptz=0) {
    double sumX = 0;
    double sumZ = 0;
    double sumXZ = 0;
    double sumZsquare = 0;  // = Sum (Z^2)
    double sumW = 0;

    for(const auto& h : cluster) {
    	if(h.flag!=PositionHit::Flag::valid) continue;
    	double errorx2=h.error.x*h.error.x;
    	double hiz=h.position.z-interceptz;
		sumX += h.position.x/errorx2;
		sumZ += hiz/errorx2;
		sumXZ += h.position.x*hiz/errorx2;
		sumZsquare += hiz*hiz/errorx2;
		sumW+=1/errorx2;
    }

    double denominator=(sumZ * sumZ - sumW * sumZsquare);
    if(std::fabs(denominator)<1E-20){
    	std::cerr<<"regressionXZ error: (sumZ * sumZ - ntot * sumZsquare)<1E-20"<<std::endl;
    	std::cerr<<"sumZ="<<sumZ<<" sumZsquare="<<sumZsquare<<" ntot="<<sumW<<std::endl;
    	std::cerr<<cluster.size()<<" hits"<<std::endl;
//    	std::cerr<<cluster.size()<<" hits on "<<cluster.getNPlanesHit()<<" planes"<<std::endl;
    	throw "(sumZ * sumZ - ntot * sumZsquare)<1E-20";
    }

    double slope1     = (sumX * sumZ - sumW * sumXZ) / denominator;
    double intersept1 = (sumZ * sumXZ - sumZsquare * sumX) / denominator;

    double sigmaIntercept2=-sumZsquare/denominator;
    double sigmaSlopeIntercept=sumZ/denominator;
    double sigmaSlope2=-sumW/denominator;
    std::vector<double> error={ sigmaSlope2, sigmaSlopeIntercept, sigmaIntercept2 };

    return FitResult2D(slope1, intersept1, error , interceptz);

}

template<class HitContainer=HoughTransformer::HitCluster>
FitResult2D regressionYZ(const HitContainer& cluster, double interceptz=0) {
    double sumY = 0;
    double sumZ = 0;
    double sumYZ = 0;
    double sumZsquare = 0;  // = Sum (Z^2)
    double sumW = 0;

    for(auto& h : cluster) {
    	if(h.flag!=PositionHit::Flag::valid) continue;
    	double errory2=h.error.y*h.error.y;//add error!
    	double hiz=h.position.z-interceptz;
		sumY += h.position.y/errory2;
		sumZ += hiz/errory2;
		sumYZ += h.position.y*hiz/errory2;
		sumZsquare += hiz*hiz/errory2;
		sumW+=1/errory2;
    }

    double denominator=(sumZ * sumZ - sumW * sumZsquare);
    if(std::fabs(denominator)<1E-20){
    	std::cerr<<"regressionYZ error: (sumZ * sumZ - ntot * sumZsquare)<1E-20"<<std::endl;
    	std::cerr<<"sumZ="<<sumZ<<" sumZsquare="<<sumZsquare<<" ntot="<<sumW<<std::endl;
    	std::cerr<<cluster.size()<<" hits"<<std::endl;
//    	std::cerr<<cluster.size()<<" hits on "<<cluster.getNPlanesHit()<<" planes"<<std::endl;
    	throw "(sumZ * sumZ - ntot * sumZsquare)<1E-20";
    }

    double slope1     = (sumY * sumZ - sumW * sumYZ) / denominator;
    double intersept1 = (sumZ * sumYZ - sumZsquare * sumY) / denominator;

    double sigmaIntercept2=-sumZsquare/denominator;
    double sigmaSlopeIntercept=sumZ/denominator;
    double sigmaSlope2=-sumW/denominator;
    std::vector<double> error={ sigmaSlope2, sigmaSlopeIntercept, sigmaIntercept2 };

    return FitResult2D(slope1, intersept1, error, interceptz );

}


template<class HitContainer>
FitResult3D regressionFit3d(const HitContainer& cluster, double interceptz) {
	return FitResult3D {
		regressionXZ(cluster, interceptz),
		regressionYZ(cluster, interceptz)
	};
}

PositionHit& calculateResidual(PositionHit& h, const FitResult3D& fit ) {
	h.residual.x=h.position.x - fit.xAt(h.position.z);
	h.residual.y=h.position.y - fit.yAt(h.position.z);
	h.residual.z=0;
	return h;
}

template<class hitContainer>
hitContainer& calculateResiduals( hitContainer& cluster, const FitResult3D& fit)  {
	for(auto& h : cluster) {
		 calculateResidual(h, fit);
	}
	return cluster;
}

//function invalidates residuals reference to hits!
HoughTransformer::HitCluster&  cutOnResiduals( HoughTransformer::HitCluster& cluster, double maxResidual ) {
	int nremoved=cluster.size();
	//no erase necessary because list version of remove_if
	cluster.remove_if( [&maxResidual](const PositionHit& h){
		bool cut= h.residual.x*h.residual.x + h.residual.y*h.residual.y > maxResidual*maxResidual;
		return cut;
	} );
	nremoved-=cluster.size();
//	std::cout<<"Removed "<<nremoved<<" hits"<<std::endl;
	return cluster;
}


TVector3 averageResidual(const HoughTransformer::HitCluster& residuals) {
	double x=0, y=0, z=0;
	for(auto& h: residuals) {
		x+=h.residual.x;
		y+=h.residual.y;
		z+=h.residual.z;
	}
	return { x/=residuals.size(),
			 y/=residuals.size(),
			 z/=residuals.size() };
}

/* TODO: remove residual class from these functions
HoughTransformer::HitCluster& cutOnResidualPulls(
		HoughTransformer::HitCluster& cluster,
		const std::vector<Residual>& residuals, double maxPullx, double maxPully) {
	auto res=residuals.begin();
	int nremoved=cluster.size();
	for(auto& h : cluster) {
		if( res->x*res->x/h.error2x > maxPullx*maxPullx) h.flag=-1;
		if( res->y*res->y/h.error2y > maxPully*maxPully) { h.flag=-2;}
		++res;
	}
	nremoved-=cluster.size();
	return cluster;
}


HoughTransformer::HitCluster& cutOnResidualPullsWithFitError(
		HoughTransformer::HitCluster& cluster,
		const FitResult3D& fit,
		double maxPullx,
		double maxPully ) {
	int nremoved=cluster.size();
	for(auto& h : cluster) {
		auto res=calculateResidual(h,fit);
		double xerror2=h.error2x+fit.XZ.error2At(h.z);
		double yerror2=h.error2y+fit.YZ.error2At(h.z);
		if( res.x*res.x/xerror2 > maxPullx*maxPullx) h.flag=-1;
		if( res.y*res.y/yerror2 > maxPully*maxPully) { h.flag=-2;}
	}
	nremoved-=cluster.size();
	return cluster;
}
*/
