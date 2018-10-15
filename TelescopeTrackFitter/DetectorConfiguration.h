/*
 * detectorConfiguration.h
 *
 *  Created on: Jun 23, 2017
 *      Author: cligtenb
 */

#ifndef DETECTORCONFIGURATION_H_
#define DETECTORCONFIGURATION_H_

#include <utility>

struct DetectorConfiguration {
	virtual ~DetectorConfiguration() {};

	virtual double xmax() const =0;
	virtual double xmin() const { return 0.; }
	virtual double ymax() const  =0;
	virtual double ymin() const { return 0.; }
	virtual double zmax() const  =0;
	virtual double zmin() const { return 0; }
	virtual double zmean() const { return (zmin()+zmax())/2; }

};

struct SimpleDetectorConfiguration : DetectorConfiguration {
	SimpleDetectorConfiguration(double minx, double maxx, double miny, double maxy, double minz, double maxz) : minx(minx), maxx(maxx), miny(miny), maxy(maxy), minz(minz), maxz(maxz) {};
	double minx, maxx, miny, maxy, minz, maxz;
	virtual ~SimpleDetectorConfiguration() {};
	virtual double xmax() const { return maxx; }
	virtual double xmin() const { return minx; }
	virtual double ymax() const { return maxy; }
	virtual double ymin() const { return miny; }
	virtual double zmax() const { return maxz; }
	virtual double zmin() const { return minz; }
};

struct TelescopeConfiguration : DetectorConfiguration {
	virtual ~TelescopeConfiguration() {}
	TelescopeConfiguration (int nPlanes, std::vector<double> planePosition, double pixelsize, int pixelColumns, int pixelRows) :
			nPlanes(nPlanes),
			planePosition(planePosition),
			pixelsize(pixelsize),
			pixelColumns(pixelColumns),
			pixelRows(pixelRows)
		{};

	int nPlanes;
	std::vector<double> planePosition; //mm
	double pixelsize;//mm
	int pixelColumns, pixelRows;

	double planexmax() const { return pixelColumns*pixelsize; }
	double planeymax() const { return pixelRows*pixelsize; }

	virtual double xmax() const { return planexmax(); }
	virtual double xmin() const { return 0.; }
	virtual double ymax() const { return planeymax(); }
	virtual double ymin() const { return 0.; }
	virtual double zmax() const { return planePosition.back(); }
	virtual double zmin() const { return planePosition.front(); }

	virtual std::pair<double,double> getCentre() const { return { (xmax()-xmin())/2, (ymax()-ymin())/2 }; };

};


const TelescopeConfiguration mimosa= {
	6, //planes
//	{0, 18.6, 37.4, 116.7, 151.1, 188.4}, //plane position from Wolf thesis
//	{0, 15.8, 31.8, 143.1, 161.55, 179.91 }, //plane positions as measured
	{0, 19.7, 38.9, 396.6, 415.8, 435.4 }, //plane positions as measured in oktober18 quad
	0.0184, //pixelsize in mm
	1152, 576 //row, column
};


#endif /* DETECTORCONFIGURATION_H_ */
