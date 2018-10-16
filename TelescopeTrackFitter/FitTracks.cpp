#include <string>
#include <iostream>

#include "TROOT.h"
#include "TGraph.h"
#include "TH1.h"

#include "/user/cligtenb/rootmacros/AllCombiner.h"

#include "DetectorConfiguration.h"
#include "TelescopeTrackFitter.h"

//#include "/user/cligtenb/rootmacros/AllCombiner.h"

#if 1 //root?
#include "TelescopeTrackFitter.cpp"
//#include "linearRegressionFit.cpp"
//#include "makeNoisyPixelMask.cpp"
#include "ResidualHistogrammer.cpp"
#endif

using namespace std;


//std::vector<std::pair<double, double>>& fixFirstPosition(std::vector<std::pair<double, double>>& means) {
//	for(int i=1; i<means.size(); ++i) {
//		means[i].first+=means[0].first;
//		means[i].second+=means[0].second;
//	}
//	means[0]=std::make_pair(0.,0.);
//	return means;
//}

void testFitTracks(std::string inputfile, std::string alignmentFile="") {

	TelescopeTrackFitter telescopeFitter(inputfile, mimosa);
	telescopeFitter.makeMask(5e3);

	telescopeFitter.fitTracks("test.root");

}

void FitTracks (std::string inputfile, int nRepeatFit=6) {

	TelescopeTrackFitter telescopeFitter(inputfile, mimosa);

	telescopeFitter.makeMask(5e3);
//	telescopeFitter.setShifts( { {0,0}, {0,0}, {0,0}, {0,0}, {0,1}, {0,0} } );

	//initialise alignment parameters
	double recursion[nRepeatFit],
		shiftx[mimosa.nPlanes][nRepeatFit],
		shifty[mimosa.nPlanes][nRepeatFit],
		rotation[mimosa.nPlanes][nRepeatFit];

	const int firstRefPlane=1, secondRefPlane=4;

	for(int i=0; i<nRepeatFit; i++) {
				cout<<"fitting "<<i<<endl;
//				if(i==4) telescopeFitter.displayEvent=true;

				if(i>=3) telescopeFitter.makeTrackHistograms=true;

				if(i<=2) telescopeFitter.selectHitForRefit=[](const PositionHit& h) {return h.chip==firstRefPlane || h.chip==secondRefPlane ;};
				else telescopeFitter.selectHitForRefit=[](const PositionHit& h) {return true;};

				if(i==1) telescopeFitter.constructLineParallelToZ=true;
				else telescopeFitter.constructLineParallelToZ=false;

				//fit tracks!
				telescopeFitter.fitTracks("residualHistograms"+to_string(i)+".root");

				auto means=telescopeFitter.getMeanResiduals();
				auto rotations=telescopeFitter.getRotations();

				means[firstRefPlane]=means[secondRefPlane]={0,0};
				rotations[firstRefPlane]=0;

				std::vector<double> case1Angle(6,0.);
				switch(i) {
				case 0:
					telescopeFitter.addToShifts( means );
					telescopeFitter.setAngles({0,0,0,0,0,0});
					break;
				case 1:
					case1Angle[secondRefPlane]=rotations[secondRefPlane];
					telescopeFitter.addToAngles( case1Angle );
					break;
				case 2:
					telescopeFitter.addToAngles( rotations );
					break;
				case 3:
					telescopeFitter.setSlopes( telescopeFitter.getSlopes() );
					telescopeFitter.maxResidual=0.05;
					break;
				default:
					telescopeFitter.addToShifts( means );
					telescopeFitter.addToAngles( rotations );
					break;
				}

				//store alignment parameters
				recursion[i]=i;
				for(int plane=0; plane<mimosa.nPlanes; ++plane) {
					shiftx[plane][i]=telescopeFitter.getShifts().at(plane).first;
					shifty[plane][i]=telescopeFitter.getShifts().at(plane).second;
					rotation[plane][i]=telescopeFitter.getAngles().at(plane);
				}
	}

	if(!telescopeFitter.displayEvent) {
		telescopeFitter.saveAlignment("align.dat");
	}

	//create and combine graphs
	std::vector<TGraph*> shiftxGraph,shiftyGraph, rotationGraph;
	for(int plane=0;plane<mimosa.nPlanes; ++plane) {
		shiftxGraph.push_back( new TGraph(nRepeatFit, recursion, shiftx[plane]) );
		shiftyGraph.push_back( new TGraph(nRepeatFit, recursion, shifty[plane]) );
		rotationGraph.push_back( new TGraph(nRepeatFit, recursion, rotation[plane]) );

		shiftxGraph.back()->SetTitle( ("plane "+to_string(plane+1)+";Iteration;Correction [mm]").c_str() );
		shiftyGraph.back()->SetTitle( ("plane "+to_string(plane+1)+";Iteration;Correction [mm]").c_str() );
		rotationGraph.back()->SetTitle( ("plane "+to_string(plane+1)+";Iteration;Correction [rad.]").c_str() );
	}

	AllCombiner<TGraph>  Xcombination("shiftsxCombined", shiftxGraph);
	AllCombiner<TGraph>  Ycombination("shiftsyCombined", shiftyGraph);
	AllCombiner<TGraph>  RotCombination("rotationCombined", rotationGraph);

	std::vector<TCanvas*> canv;
	for(AllCombiner<TGraph>* comb : {&Xcombination, &Ycombination, &RotCombination}) {
		comb->setStyle(7);
		canv.push_back( comb->createCombined() );
	}

}


//redirect to root main function
int main(int argc, const char* argv[]) {
	gROOT->ProcessLine(".L Hit.h+"); //compile hit
	if(argc>1) FitTracks(argv[1]);
}
