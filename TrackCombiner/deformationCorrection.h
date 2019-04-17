/*
 * deformationCorrection.h
 *
 *  Created on: Nov 30, 2018
 *      Author: cligtenb
 */

#ifndef TRACKCOMBINER_DEFORMATIONCORRECTION_H_
#define TRACKCOMBINER_DEFORMATIONCORRECTION_H_

#include <vector>


#include "TF1.h"
#include "TF2.h"

//peters correction
//double deformationCorrection(int ichip, double x) {
//// fit parameters for chip 0 till 3
// std::vector <double>  p {11, 1.52397, 1.2152, -0.686208, 1, 1.39945, -0.00179486,
//                          11.5, 0.443004, 0.545847, -0.671393, 0.233943, 2.83167, 0.0402569,
//                          12.5, -1.01292, 0.81062, 0.564847, -0.597693, 2.31683, 0.025552,
//                          12, -9.97675, 1.92863, 0.0513599, -9.91428, 0.979708, 0.00343746 };
// if(ichip>4) return 0;
// if(ichip<0) return 0;
// int i = 7*ichip;
// return -p[i+1]/(1+(x-p[i])*(x-p[i])/p[i+2]/p[i+2])+p[i+4]/(1+(x-p[i]-p[i+3])*(x-p[i]-p[i+3])/p[i+2]/p[i+2]/p[i+5])+p[i+6];
//}


//fitted 2x breit-wigners
//double deformationCorrection(int ichip, double x) {
//	static auto correction=new TF1("correction", "-[p1]/(1+std::pow((x-[p0])/[p2], 2)) + [p4]/(1+std::pow((x-[p3])/[p5],2))+[p6]");
//
//	const std::vector<std::vector<double>> p{
//		{13.34, -2.01364, 1.18425, 13.6036, -2.40334, 1.0573, -0.0177269},
//		{13.032, -0.661975, 1.60116, 16.128, -3.64937, 1.5225, 0.0927064},
//		{14.6021, -4.666, 0.921731, 14.6671, -4.2054, 1.03022, 0.0100087},
//		{14.9556, -1.75714, 0.703961, 15.03, -1.5939, 0.851984, 0.00892923} };
//
//	for(int i=0;i<7;i++) correction->SetParameter(i,p[ichip][i]);
//	return correction->Eval(x);
//}


//static auto correction=new TF1("correction",
//		"breitwigner(0)+breitwigner(3)+breitwigner(6)+breitwigner(9)+[offset]");

const std::vector<std::vector<double>> deformationCorrectionParameters{
	{2.61443, 0.100941, 1.61043, -1.52508, 0.926549, 1.86461, 1.36493, 13.126, 1.96334, -3.15446, 14.3033, 1.87523, -0.000113656},
	{2.2285, -0.0143355, 1.65785, -1.16026, 0.995182, 1.79715, 2.64069, 13.2987, 2.55817, -6.25284, 14.6259, 2.36398, 0.000223841},
	{4.19523, 14.0504, 2.02393, -2.35141, 14.9244, 2.3167, 1.72633, 27.3346, 1.94821, -3.01411, 28.2367, 1.71122, 0.00287906},
	{6.24837, 13.7489, 2.48634, -2.89738, 15.1342, 2.68989, 1.65655, 27.4766, 2.02193, -3.22435, 28.5387, 1.84855, 0.00445736}};


double deformationCorrection(int ichip, double x, const std::vector<std::vector<double>>& p=deformationCorrectionParameters) {
	static auto correction=new TF1("correction",
			"breitwigner(0)+breitwigner(3)+breitwigner(6)+breitwigner(9)+[offset]");

	for(std::size_t i=0;i<p[ichip].size();i++) correction->SetParameter(i,p[ichip][i]);

	return correction->Eval(x);
}

double deformationCorrection2D(int ichip, double x, double y, const std::vector<std::vector<double>>& p=deformationCorrectionParameters) {
	auto correction=new TF2("correction",
			"(1+[b0]*y+[c0]*y*y+[d0]*y*y*y+[e0]*y*y*y*y)*breitwigner(0)+(1+[b1]*y+[c1]*y*y+[d1]*y*y*y+[e1]*y*y*y*y)*breitwigner(3) "
			"+"
			"(1+[b2]*y+[c2]*y*y+[d2]*y*y*y+[e2]*y*y*y*y)*breitwigner(6)+(1+[b3]*y+[c3]*y*y+[d3]*y*y*y+[e3]*y*y*y*y)*breitwigner(9)+[offset]");

	for(std::size_t i=0;i<p[ichip].size();i++) correction->SetParameter(i,p[ichip][i]);

	return correction->Eval(x,y);
}

//double deformationCorrection(int ichip, double x) {
//static auto correction=new TF1("correction", "[p0]+x*[p1]+x*x*[p2]+x*x*x*[p3]");
//
//	const std::vector<std::vector<double>> p{
//		{-0.226695, 0.0932323, -0.0136263, 0.000614136},
//		{-0.211077, 0.0851089, -0.0125936, 0.000576097},
//		{-4.51686, 0.656857, -0.0316985, 0.000507013},
//		{-5.15714, 0.735173, -0.0348212, 0.000547832}
//	};
//
//	for(int i=0;i<correction->GetNpar();i++) correction->SetParameter(i,p.at(ichip).at(i));
//
//	struct { double min, max; } chipRange[4] = { {1.5,13}, {1.5,13}, {15.5,27}, {15.5,27} };
//	x=std::min(x, chipRange[ichip].max);
//	x=std::max(x, chipRange[ichip].min);
//
//	return correction->Eval(x);
//}

#endif /* TRACKCOMBINER_DEFORMATIONCORRECTION_H_ */
