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


static auto correction=new TF1("correction",
		"breitwigner(0)+breitwigner(3)+breitwigner(6)+breitwigner(9)+[offset]");

double deformationCorrection(int ichip, double x) {
	static auto correction=new TF1("correction",
			"breitwigner(0)+breitwigner(3)+breitwigner(6)+breitwigner(9)+[offset]");

	const std::vector<std::vector<double>> p{
		{7.09713, 0.430926, 1.79566, -6.22156, 0.684718, 1.89835, 1.84525, 13.2064, 2.07817, -4.34259, 14.382, 2.12304, 0.0103622},
		{2.48729, -0.229976, 1.33571, -0.911999, 1.19903, 1.69983, 9.1324, 13.5909, 2.26925, -11.5126, 13.9285, 2.25019, 0.00336787},
		{13.7212, 13.5094, 3.8325, -7.44684, 14.82, 3.62685, 5.96425, 27.6485, 2.07897, -7.02432, 27.9569, 1.94964, -0.033238},
		{13.7446, 13.8724, 2.85009, -10.3855, 14.4876, 3.33408, 3.03169, 27.4861, 2.27521, -5.85349, 28.5317, 2.54876, 0.0326382}};

	for(std::size_t i=0;i<p[ichip].size();i++) correction->SetParameter(i,p[ichip][i]);

	return correction->Eval(x);
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
