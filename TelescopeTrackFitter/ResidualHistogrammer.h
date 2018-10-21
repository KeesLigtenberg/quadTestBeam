/*
 * ResidualHistogrammer.h
 *
 *  Created on: Jun 15, 2017
 *      Author: cligtenb
 */

#ifndef RESIDUALHISTOGRAMMER_H_
#define RESIDUALHISTOGRAMMER_H_
#include <vector>
#include <string>

#include "TFile.h"
#include "TH1.h"
#include "TProfile.h"
#include "TProfile2D.h"

#include "linearRegressionFit.h"
#include "DetectorConfiguration.h"

#include "getMeanFromGausFit.h"

class ResidualHistogrammer {
public:
	ResidualHistogrammer(std::string outputFileName, const TelescopeConfiguration& detector);
	virtual ~ResidualHistogrammer();

	void fill(const PositionHit&, const std::pair<double, double>& rotationPoint);
	void fill(const HoughTransformer::HitCluster&, const std::vector<std::pair<double, double> >& rotationPoints );
	void fill(const HoughTransformer::HitCluster& residuals, const std::pair<double, double>& rotationPoint);
	void fill(const HoughTransformer::HitCluster& residualVector) {
		std::vector< std::pair<double,double> > rotationPoints(detector.nPlanes, detector.getCentre());
		fill(residualVector, rotationPoints );
	}


	std::vector< std::pair<double, double> > getMeansOfPlanes();
	std::vector< double > getRotationOfPlanes();

	TFile outputFile;
	const TelescopeConfiguration& detector;

	struct PlaneHistograms {
		static int n;
		PlaneHistograms(const TelescopeConfiguration& det) :
			xResidual( ("xResidual_"+std::to_string(n)).c_str(), ("xResidual_"+std::to_string(n)).c_str(), 40, -0.2, 0.2),
			yResidual( ("yResidual_"+std::to_string(n)).c_str(), ("yResidual_"+std::to_string(n)).c_str(), 40, -0.2, 0.2),
			zRotation( ("zRotation_"+std::to_string(n)).c_str(), ("zRotation_"+std::to_string(n)).c_str(), 40, -0.05, 0.05),
			xResidualByPixel( ("xResidualByPixel_"+std::to_string(n)).c_str(), ("xResidualByPixel_"+std::to_string(n)).c_str(), 20, det.xmin(), det.xmax(), 20, 0, det.ymax() ),
			yResidualByPixel( ("yResidualByPixel_"+std::to_string(n)).c_str(), ("yResidualByPixel_"+std::to_string(n)).c_str(), 20, det.xmin(), det.xmax(), 20, 0, det.ymax() ),
//			xRotation( ("xRotation_"+std::to_string(n)).c_str(), ("xRotation_"+std::to_string(n)).c_str(), 40, -0.05, 0.05),
//			yRotation( ("yRotation_"+std::to_string(n)).c_str(), ("yRotation_"+std::to_string(n)).c_str(), 40, -0.05, 0.05),
			zRotationByPixel( ("zRotationByPixel_"+std::to_string(n)).c_str(), ("zRotationByPixel_"+std::to_string(n)).c_str(), 20, det.xmin(), det.xmax(), 20, 0, det.ymax() ),
			hitmap( ("hitmap_"+std::to_string(n)).c_str(), ("hitmap_"+std::to_string(n)).c_str(), det.pixelColumns, 0, det.pixelColumns, det.pixelRows, 0, det.pixelRows )
		{ ++n; }
		std::pair<double,double> getMeansFromFit();
		double getRotationFromFit();
		TH1D xResidual, yResidual;
		TH1D zRotation;// xRotation, yRotation;
		TH2D hitmap;
		TProfile2D xResidualByPixel, yResidualByPixel, zRotationByPixel;
	};
	std::vector<PlaneHistograms> planeHist;
	TProfile xResidualByToT;
	TH1D ToT;

};


class TrackHistogrammer {
public:
	TrackHistogrammer(const DetectorConfiguration& detector);
	void fill(const FitResult3D& entry);
private:
	//	TH1D phi, d0, tanLambda, z0;
	TH1D slope1, slope2, intercept1, intercept2;
	const DetectorConfiguration& detector;
};


#endif /* RESIDUALHISTOGRAMMER_H_ */
