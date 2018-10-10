//root macro for quick plots

#include <cstdlib>
#include <iostream>

#include "TSystem.h"

#include "/user/cligtenb/rootmacros/getObjectFromFile.h"
#include "/user/cligtenb/rootmacros/getHistFromTree.h"
#include "/user/cligtenb/rootmacros/AllCombiner.h"

using namespace std;


const int nChips=4;
struct shift {double x,y;};
std::array<shift,nChips> chipShifts={{ {14.014,0.012}, {-0.026, 10.659}, {14.237, 10.659}, {28.273, 0.006} }};

void animateHitmap(std::string filename="tree.root", std::string treename="data") {
	auto tree=getObjectFromFile<TTree>(treename, filename);
	
	TH2D* axis=new TH2D("axis", ";x [mm];y [mm]", 1,0,28,1,-14,25);
	new TCanvas("hitMapCanvas", "Canvas for hit map", 28*20,38*20);
	gStyle->SetOptStat(0);

	int nEntries=tree->GetEntries();
	const int entriesPerFrame=1;
	for(int i=0; i<nEntries; i+=entriesPerFrame) {
		axis->Draw();
		tree->SetMarkerStyle(7);	
		const int nChips=4;
		for(int j=0; j<nChips; j++) {
			
			const std::array<int,nChips> colorNumbers{ kRed-8, kBlue-8, kGreen-8, kMagenta-8};

			tree->SetMarkerColor(colorNumbers[j]);
			auto js=to_string(j);
			auto dir=to_string( (j==0 || j==3) ? -1 : 1 );
			tree->Draw( (dir+"*chip"+js+".row*.055+"+to_string(chipShifts[j].y)+":"+dir+"*chip"+js+".column*.055+"+to_string(chipShifts[j].x)).c_str(), "", "same", entriesPerFrame, i);
		}
		cout<<i<<"\n";
		gPad->Update();
		gPad->SetTicks(1,1);
//		gSystem->Sleep(50);
		if(cin.get()=='q') break;
	}

}


		
void drawHitmap(std::string filename="tree.root", std::string treename="data") {
	auto tree=getObjectFromFile<TTree>(treename, filename);
	
	TH2D* axis=new TH2D("axis", ";x [mm];y [mm]", 1000,0,28,1000,-14,25);
	new TCanvas("hitMapCanvas", "Canvas for hit map", 28*20,38*20);
	gStyle->SetOptStat(0);

	int nEntries=tree->GetEntries();
	axis->Draw();
	tree->SetMarkerStyle(7);	
	for(int j=0; j<nChips; j++) {
		auto dir=to_string( (j==0 || j==3) ? -1 : 1 );
		const std::array<int,nChips> colorNumbers{ kRed-8, kBlue-8, kGreen-8, kMagenta-8};
		tree->SetMarkerColor(colorNumbers[j]);
		auto js=to_string(j);
		tree->Draw( (dir+"*chip"+js+".row*.055+"+to_string(chipShifts[j].y)+":"+dir+"*chip"+js+".column*.055+"+to_string(chipShifts[j].x)+">>+axis").c_str(), "", "colzsame");

	}
	gPad->Update();
	gPad->SetTicks(1,1);
}
/*

void plotSpot( std::string filename="tree.root", std::string treename="data") {
	auto tree=getObjectFromFile<TTree>(treename, filename);
	auto canv=new TCanvas("spot", "Measured-Expected", 1200,300);
	canv->Divide(nChips,1);
	for(int i=0;i<nChips;i++) {
		auto dir=to_string( (i==0 || i==3) ? -1 : 1 );
		
		auto spot=getHistFromTree(
			*tree, 
			dir+"*row*0.055-laser.y+"+to_string(chipShifts[i].y) 
			+":"+
			dir+"*-column*0.055-laser.x+"+to_string(chipShifts[i].x) 
		 ,"chip=="+to_string(i), "colzgoff" );
		spot->SetTitle(";x_{Chip}-x_{Laser} [mm];y_{Chip}-y_{Laser} [mm]");
		canv->cd(i+1);
		spot->Draw();
		gPad->SetTicks(1,1);
		increaseAxisSize(spot);
	}
}

void plotSpotPixels( std::string filename="tree.root", std::string treename="data") {
	auto tree=getObjectFromFile<TTree>(treename, filename);
	auto canv=new TCanvas("spot", "Measured-Expected", 1000,1000);
	canv->Divide(2,2);
	for(int i=0;i<nChips;i++) {
		auto dir=to_string( (i==0 || i==3) ? -1 : 1 );
		
		auto spot=getHistFromTree(
			*tree, 
			dir+"*row+(-laser.y+"+to_string(chipShifts[i].y)+")/.055" 
			+":"+
			dir+"*-column+(-laser.x+"+to_string(chipShifts[i].x)+")/.055"
			,"chip=="+to_string(i),
			"chip"+to_string(i)+"(100,-50,50,100,-50,50)",
			"colzgoff" );
		spot->SetTitle(";x_{Chip}-x_{Laser} [columns];y_{Chip}-y_{Laser} [rows]");
		canv->cd(i+1);
		spot->Draw("colz");
		gPad->SetTicks(1,1);
		increaseAxisSize(spot);
	}
}

void plotDriftTime( std::string filename="tree.root", std::string treename="data") {
	auto tree=getObjectFromFile<TTree>(treename, filename);
	auto drift=getHistFromTree(*tree, "driftTime/4096*25E-3");
	drift->SetTitle(";drift time [#mus];Hits");
	drift->Draw();
	gPad->SetTicks(1,1);
	increaseAxisSize(drift);
}

void plotToT( std::string filename="tree.root", std::string treename="data") {
	auto tree=getObjectFromFile<TTree>(treename, filename);
	auto ToT=getHistFromTree(*tree, "ToT*25E-3", "", "ToT(60,0,1.5)", "");
	ToT->SetTitle(";ToT [#mus];Hits");
	ToT->Draw();
	gPad->SetTicks(1,1);
	increaseAxisSize(ToT);
}


void plotNHits( std::string filename="tree.root", std::string treename="data") {
	auto tree=getObjectFromFile<TTree>(treename, filename);
	auto ToT=getHistFromTree(*tree, "Length$(hits)");
	//auto ToT=getHistFromTree(*tree, "Length$(hits)", "(chip==3 && fabs(0.055*row+fY-15.5)<1.5 && fabs(0.055*column-fX+10)<1.5)*1./Length$(hits)"); does not work yet
	ToT->SetTitle(";Number of hits per laser pulse;");
	ToT->Draw("HIST");
	gPad->SetTicks(1,1);
	increaseAxisSize(ToT);
}
*/


