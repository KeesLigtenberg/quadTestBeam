/*
 * deformationCorrection.h
 *
 *  Created on: Nov 30, 2018
 *      Author: cligtenb
 */

#ifndef TRACKCOMBINER_DEFORMATIONCORRECTION_H_
#define TRACKCOMBINER_DEFORMATIONCORRECTION_H_

#include <vector>


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

double deformationCorrection(int ichip, double x) {
	static auto correction=new TF1("correction", "-[p1]/(1+std::pow((x-[p0])/[p2], 2)) + [p4]/(1+std::pow((x-[p3])/[p5],2))+[p6]");

	const std::vector<std::vector<double>> p{
		{13.34, -2.01364, 1.18425, 13.6036, -2.40334, 1.0573, -0.0177269},
		{13.032, -0.661975, 1.60116, 16.128, -3.64937, 1.5225, 0.0927064},
		{14.6021, -4.666, 0.921731, 14.6671, -4.2054, 1.03022, 0.0100087},
		{14.9556, -1.75714, 0.703961, 15.03, -1.5939, 0.851984, 0.00892923} };

	for(int i=0;i<7;i++) correction->SetParameter(i,p[ichip][i]);
	return correction->Eval(x);
}

#endif /* TRACKCOMBINER_DEFORMATIONCORRECTION_H_ */
