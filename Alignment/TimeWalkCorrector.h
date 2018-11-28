/*
 * TimeWalkCorrector.h
 *
 *  Created on: Jul 20, 2018
 *      Author: cligtenb
 */

#ifndef LASERDATAFITTER_TIMEWALKCORRECTOR_H_
#define LASERDATAFITTER_TIMEWALKCORRECTOR_H_

#include <memory>

#include "TF1.h"
#include "TTree.h"
#include "TH2D.h"
#include "TROOT.h"
#include "TPad.h"

#include "../TrackFitter/PositionHit.h"
#include "AlignmentHolder.h"

struct TimeWalkCorrector : AlignmentHolder {
public:
	TimeWalkCorrector() : AlignmentHolder("TIMEWALKPARAMETERS") {}

	std::vector<PositionHit>&  correct(std::vector<PositionHit>& spaceHit) const;
	PositionHit&  correct(PositionHit& spaceHit) const;
	double getCorrection(double ToT) const;

	void update(TH2D* zResidualByToT);

	void writeParameters(std::ostream&);
	void readParameters(std::istream&);

	const std::vector<double>& getParameters() const {return params;}

private:
	double minToT{0.15};
	std::unique_ptr<TF1> fun{};
	std::vector<double> params{};
	std::string funString{};
};


double TimeWalkCorrector::getCorrection(double ToT) const {
	ToT*=0.025;
	if(ToT<minToT) ToT=minToT;
	if(fun)	return fun->Eval(ToT);
	else {
		std::cerr<<"error: function was not defined!\n";
		return 0;
	}
}

std::vector<PositionHit>& TimeWalkCorrector::correct(std::vector<PositionHit>& spaceHit) const {
//	spaceHit.erase(
//			std::remove_if(spaceHit.begin(), spaceHit.end(), [this](const PositionHit&h) {return (h.ToT*0.025) <minToT;} ),
//			spaceHit.end()
//	);
	for(auto& h : spaceHit) {
		h=correct(h);
	}
	return spaceHit;
}
PositionHit&  TimeWalkCorrector::correct(PositionHit& h) const {
	if(h.ToT/40.<minToT) {h.flag=PositionHit::Flag::lowToT;};
	h.position.z=h.position.z-getCorrection( h.ToT );
	return h;
}

void TimeWalkCorrector::writeParameters(std::ostream& output) {
	output<<minToT<<" "<<params.size()<<"\n";
	output<<funString<<"\n";
	for(const auto& x : params) {
		output<<x<<" ";
	}
	output<<"\n";
}

namespace {
	void checkStream( std::istream& input , std::string message="") {
		if(not input.good() ) {
			std::cerr<<"Failed to read parameters <<" << message <<"\n";
			throw 1;
		}
	}
}
void TimeWalkCorrector::readParameters(std::istream& input) {
	int n;
	checkStream(input, "TimeWalkCorrector: header");
	input>>minToT>>n;
	checkStream(input, "TimeWalkCorrector: minToT n");
	if(input.peek()=='\n') input.get();
	getline(input,funString);
	checkStream(input, "TimeWalkCorrector: funstring");
	fun=std::unique_ptr<TF1>( new TF1("timewalkCorrection", funString.c_str()) );
	for(int i=0; i<n; i++) {
		double p;
		input>>p;
		params.push_back(p);
		checkStream(input, "TimeWalkCorrector: parameter "+std::to_string(i) );
		fun->SetParameter(i, p);
	}
	//first parameter should be offset, so set to zero
	fun->SetParameter(0,0);
}

void TimeWalkCorrector::update(TH2D* th2) {
	TF1* gausRangeTW=new TF1("gausRangeTW","gaus(0)", -5,5);
	gausRangeTW->SetParameters(4E4,0.05,0.22);
	th2->FitSlicesY(gausRangeTW,0,-1,5);
	std::string meanHistName=th2->GetName()+std::string("_1");
	auto means=dynamic_cast<TH1*>( gDirectory->Get(meanHistName.c_str()) );
	if(!means) { std::cerr<<"failed to retrieve result from fitslicesy()\n"; return; };

//	means->Fit(fun.get(), "QS", "", 0,2.5);
	means->Fit(fun.get(), "QS", "", minToT, 2.5);
	if(fun->GetNpar()!= int(params.size())) {std::cerr<<"number of parameters in fit and in file does not match!\n"; return; };
	params[0]=1; //first parameter is offset, so set to zero!
	for(int i=1; i<fun->GetNpar(); i++) {
		std::cout<<"update parameter "<<fun->GetParName(i)<<" from "<<params[i]<<" to "<<fun->GetParameter(i)<<"\n";
		params[i]=fun->GetParameter(i);
	}
	gPad->Update(); std::cin.get();
}


#endif /* LASERDATAFITTER_TIMEWALKCORRECTOR_H_ */
