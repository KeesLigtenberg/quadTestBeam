//root macro
#include <cstdlib>
#include <vector>

#include "TPolyLine.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TH3.h"
#include "TF2.h"

#include "/user/cligtenb/rootmacros/getObjectFromFile.h"
#include "/user/cligtenb/rootmacros/getHistFromTree.h"
#include "/user/cligtenb/rootmacros/AllCombiner.h"
#include "/user/cligtenb/rootmacros/histogramOperations.h"
#include "/user/cligtenb/rootmacros/StatsWrapper.h"
#include "/user/cligtenb/rootmacros/CombineHistogramsFromTree.h"
#include "/user/cligtenb/rootmacros/CombineHistogramsFromFiles.h"

//#include "../rootmacros/getObjectFromFile.h"
//#include "../rootmacros/getHistFromTree.h"
//#include "../rootmacros/AllCombiner.h"
//#include "../rootmacros/histogramOperations.h"
//#include "../rootmacros/StatsWrapper.h"


#include "Alignment/Alignment.h"
#include "TelescopeTrackFitter/linearRegressionFit.h"
#include "TrackCombiner/deformationCorrection.h"

using namespace std;

const int nChips=4;
std::vector<std::string> chipDirectories={"chip0", "chip1","chip2","chip3"};
std::vector<int> chipMapping={2,1,3,0};
std::vector<int> chipReverseMapping={4,2,1,3};//for canv->cd

void drawChipEdges(std::string alignFile="align.dat") {
	Alignment	alignment(alignFile);
	alignment.drawChipEdges();
	gStyle->SetOptStat(0);
}

void drawChipEdgesLocal(std::string alignFile="align.dat") {
	Alignment	alignment(alignFile);
	alignment.drawChipEdges(false);
	gStyle->SetOptStat(0);
}

void combineHistogramsForChips(std::string histName="nHits", std::string fileName="combinedFit.root", bool addQuad=false) {
	HistogramCombiner combination(histName+"combination");

	if(addQuad) {
		combination.setStyle(3);
		auto quadHist=getObjectFromFile<TH1>("quad/"+histName, fileName);
		if(histName=="nHits") {
			quadHist->Rebin(2);
			quadHist=convertToTH1D(quadHist);
		}
		combination.add( quadHist , "All chips" );
	}

	for(int i=0; i<4; i++) {
		auto hist=getObjectFromFile<TH1>(chipDirectories[i]+"/"+histName, fileName);
		if(histName=="nHits") {
			hist->Rebin(2);
			hist=convertToTH1D(hist);
		}
//		hist->Draw(); gPad->Update(); cin.get();
		combination.add( hist , "chip "+std::to_string(i+1) );
	}
	combination.normalise();
	combination.createCombined();
}
void combineDrawHistogramsForChips(std::string expression, std::string cut="1", std::string drawOption="", std::string histName="Hist", std::string fileName="combinedFit.root", std::string treeName="fitResults") {
	auto tree = getObjectFromFile<TTree>(treeName, fileName);
	HistogramCombiner combination(histName+"Combination");
	for(int i=0;i<4;i++) {
		auto hist=getHistFromTree(tree, expression, cut+" && chip=="+to_string(i), "chip"+to_string(i)+histName, drawOption.c_str(),1E6);
		hist->SetTitle(("Chip "+to_string(i+1)).c_str());
		combination.add(hist);
	}
	combination.createCombined();
}

void combineHistogramsForChipsOnPads(std::string histName="nHits", std::string fileName="fitted.root") {
	auto canv=new TCanvas(("canv"+histName).c_str(), "combined canvas", 1000,1000);
	canv->Divide(2,2);
	for(int i=0;i<4;i++) {
		auto hist=getObjectFromFile<TH1>(chipDirectories[i]+"/"+histName, fileName);
		hist->SetTitle(("Chip "+to_string(i+1)).c_str());
		canv->cd(i+1);
		hist->Draw("colz");
		getMeanFromGausFitAroundMean(*hist,0.8,0.2);
	}
}

//draw four histograms on 4 pads in one canvas
void combineDrawPadsForChips(std::string expression, std::string cut, std::string drawOption, std::string histName="Hist", std::string fileName="fitted.root", std::string treeName="fitResults") {
	auto tree = getObjectFromFile<TTree>(treeName, fileName);
	auto canv=new TCanvas("canvas", "combined canvas", 1000,1000);
	canv->Divide(2,2);
	for(int i=0;i<4;i++) {
		canv->cd(i+1);
		auto hist=getHistFromTree(tree, expression, cut+" && chip=="+to_string(i), "chip"+to_string(i)+histName, drawOption.c_str());
		hist->SetTitle(("Chip "+to_string(i+1)).c_str());
	}
}

//only accapted hits
void plotHitMap(std::string file="combinedFit.root", std::string alignFile="align.dat", std::string hitMapName="local/positionHitMap") {
	new TCanvas("hitmap", "hitmap", 850,1000)	;

	//hits
	auto first=true;

//	for(const auto& dir:chipDirectories)
	std::string dir="quad";
	{
		auto hist=getObjectFromFile<TH2>(dir+"/"+hitMapName, file);
		if(first) {
			hist->SetTitle(";x-position [mm];y-position [mm]; Hits");
			hist->Draw("colz");
			gPad->SetMargin(0.1,0.15,0.1,0.05);
			gPad->SetTicks(1,1);

			hist->GetYaxis()->SetTitleOffset(1.3);
			hist->GetXaxis()->SetTitleOffset(1.1);
			hist->GetZaxis()->SetTitleOffset(1.3);
			
		} else {
			hist->Draw("colsame");
		}
		hist->Draw( first ? (first=false, "colz") : "colsame" );
	}

	//chip edges
//	drawChipEdges(alignFile);
	drawChipEdgesLocal(alignFile);
}

//this is all hits
void plotHitMapFromDraw(std::string file="fitted.root", std::string alignFile="align.dat") {

	new TCanvas("hitmap", "hitmap", 850,1000)	;
	auto tree=getObjectFromFile<TTree>("fitResults",file);
	auto hist=getHistFromTree(*tree, "position.y:position.x", "", "h(564,-5,27,745,173,214)", "colz",1E6);

	hist->SetTitle(";x-position [mm];y-position [mm]; Hits");
	gPad->SetMargin(0.1,0.15,0.1,0.05);
	gPad->SetTicks(1,1);

	hist->GetYaxis()->SetTitleOffset(1.3);
	hist->GetXaxis()->SetTitleOffset(1.1);
	hist->GetZaxis()->SetTitleOffset(1.3);

	drawChipEdges(alignFile);

}

void combineHistogramsWithMean(std::string histName="xResidual", std::string fileName="combinedFit.root") {
	HistogramCombiner combination(histName+"combination");
	for(int i=0; i<4; i++) {
		auto hist=getObjectFromFile<TH1>(chipDirectories[i]+"/"+histName, fileName);
//		hist->Draw(); gPad->Update(); cin.get();
		auto meanstr=std::to_string(int(1000*hist->GetMean()));
		combination.add( hist , "chip "+std::to_string(i+1)+" (mean "+meanstr+" #mum)" );
	}
	//combination.normalise();
	combination.createCombined();
}

void drawToT(std::string fileName="combinedFit.root") {
	combineDrawHistogramsForChips("ToT*25E-3", "fabs(ToT*25E-3)<2.5", "", "hist(100,0,2.5)", fileName);
}

void combineColumnToACorrection(std::string histName="zResidualByPixel", std::string fileName="combinedFit.root") {
	HistogramCombiner combination(histName+"combination");
	for(int i=0; i<4; i++) {
		auto prof=getObjectFromFile<TProfile2D>(chipDirectories[i]+"/"+histName, fileName);
//		hist->Draw(); gPad->Update(); cin.get();
		auto hist=prof->ProfileX();
		combination.add( hist , "chip "+std::to_string(i+1) );
	}
	//combination.normalise();
	combination.createCombined();

	new TCanvas("quad","quad", 500,500);
	auto prof=getObjectFromFile<TProfile2D>("quad/"+histName, fileName);
	prof->ProfileX()->Draw();
}

void plotDeformations(std::string x="x", std::string file="fitted.root", std::string alignFile="../align.dat") {
	new TCanvas(("deformations"+x).c_str(), ("deformations"+x).c_str(), 900,1000)	;
	
	auto tree=getObjectFromFile<TTree>("fitResults", file);
	//auto deformation=getHistFromTree(tree, "hitAverage."+x+"-laser."+x+":laser.y:laser.x", "hitAverage.x>0", "deformation(29,9.5,38.5,35, 1.5, 40.5)", "profcolz");
	auto deformation=getHistFromTree(tree, "residual."+x+":laser.y:laser.x", "flag>0 && hitAverage.x>0", "deformation(29,11.5,40.5,35, 1.5, 40.5)", "profcolz");
	deformation->SetBins(29,11.5,40.5,39, 1.5, 40.5);
	setMinMax((TProfile2D*) deformation, -0.2,0.2);
	deformation->SetTitle((";x [mm];y [mm];Mean "+x+"-residual [mm]").c_str());
	gPad->SetMargin(0.1,0.2,0.1,0.05);
	gPad->SetTicks(1,1);

	deformation->GetYaxis()->SetTitleOffset(1.3);
	deformation->GetXaxis()->SetTitleOffset(1.1);
	deformation->GetZaxis()->SetTitleOffset(1.6);

	drawChipEdges(alignFile);
}
void plotDeformationsXYZ(std::string file="fitted.root", std::string alignFile="../align.dat") {
	plotDeformations("x", file,alignFile);
	plotDeformations("y", file,alignFile);
	plotDeformations("z", file,alignFile);
}
//draw deformations as pixel coordinates
void plotDeformationsPixels(std::string x="x", std::string fileName="fitted.root", std::string treeName="fitResults") {
	auto tree = getObjectFromFile<TTree>(treeName, fileName);
	auto canv=new TCanvas(("canvasDeformationPixels"+x).c_str(), "combined canvas", 1200,1000);
	canv->Divide(2,2);
	for(int i=0;i<4;i++) {
		canv->cd(i+1);
		auto chipHist=new TProfile2D(("chipHist"+to_string(i)).c_str(), "-residuals per pixel", 14,0,256, 14,0,256);
		getHistFromTree(tree, "residual."+x+":row:column", "flag>0 && hitAverage.x>0 && chip=="+to_string(i), "chipHist"+to_string(i), "profcolz0");
		//chipHist->SetBins(256,0,256, 256,0,256);
		chipHist->SetMinimum(-0.2);
		chipHist->SetMaximum(0.2);
		//setMinMax((TProfile2D*) chipHist,-0.2,0.2);
		chipHist->SetTitle(("Chip "+to_string(i+1)+";Columns;Rows;Mean "+x+"-residual [mm]").c_str());
		gPad->SetMargin(0.1,0.2,0.1,0.1);
		gPad->SetTicks(1,1);
		chipHist->GetYaxis()->SetTitleOffset(1.3);
		chipHist->GetXaxis()->SetTitleOffset(1.1);
		chipHist->GetZaxis()->SetTitleOffset(1.6);
		
		//chipHist->GetXaxis()->SetLabelSize(0.05);
		chipHist->GetXaxis()->SetNdivisions(408, false);
		//chipHist->GetYaxis()->SetLabelSize(0.05);
		chipHist->GetYaxis()->SetNdivisions(408, false);
	}
	gStyle->SetOptStat(0);
}

//draw deformations as pixel coordinates
void plotToTPerPixel(std::string fileName="fitted.root", std::string treeName="fitResults") {
	auto tree = getObjectFromFile<TTree>(treeName, fileName);
	auto canv=new TCanvas("canvas", "combined canvas", 1200,1000);
	canv->Divide(2,2);
	for(int i=0;i<4;i++) {
		canv->cd(i+1);
		auto chipHist=new TProfile2D(("chipHist"+to_string(i)).c_str(), "-residuals per pixel", 32,0,256, 32,0,256);
		getHistFromTree(tree, "ToT/40.:row:column", "flag>0 && hitAverage.x>0 && chip=="+to_string(i), "chipHist"+to_string(i), "profcolz0");
		//chipHist->SetBins(256,0,256, 256,0,256);
		chipHist->SetMinimum(0.3);
		chipHist->SetMaximum(0.7);
		//setMinMax((TProfile2D*) chipHist,0,2);
		chipHist->SetTitle(("Chip "+to_string(i+1)+";Columns;Rows;Mean ToT [#mus]").c_str());
		gPad->SetMargin(0.1,0.2,0.1,0.1);
		gPad->SetTicks(1,1);
		chipHist->GetYaxis()->SetTitleOffset(1.3);
		chipHist->GetXaxis()->SetTitleOffset(1.1);
		chipHist->GetZaxis()->SetTitleOffset(1.6);
		
		//chipHist->GetXaxis()->SetLabelSize(0.05);
		chipHist->GetXaxis()->SetNdivisions(408, false);
		//chipHist->GetYaxis()->SetLabelSize(0.05);
		chipHist->GetYaxis()->SetNdivisions(408, false);
	}
	gStyle->SetOptStat(0);
}

void plotHitAverage(std::string file="fitted.root") {
	new TCanvas("hitAverage", "hitAverage", 750,1000)	;
	TTree* tree=getObjectFromFile<TTree>("fitResults",file);

	//hit positions
	TGraph* hitGraph=getGraphFromTree(*tree,"hitAverage.y:hitAverage.x", "hitAverage.x>0 && nHitsPassed>3");
	//hitGraph->SetMarkerStyle(7);
	hitGraph->SetTitle(";x [mm];y [mm]");
	hitGraph->Draw("AP");

	//chip edges
	drawChipEdges();
}

void plotResidualsTimeWalk( std::string filename="fitted.root" ) {

	auto uncorrected=getObjectFromFile<TH2>("quad/zResidualByToT", filename);
	TH1* uncorredtedHist=uncorrected->ProjectionY();
	auto corrected=getObjectFromFile<TH2>("quad/zResidualByToTCorrected", filename);
	TH1* correctedHist=corrected->ProjectionY();
	uncorredtedHist=makeXShifted(uncorredtedHist, -uncorredtedHist->GetMean() );

	HistogramCombiner combination("zResiduals time walk");
	combination.add(uncorredtedHist, "z-residuals without time walk correction;z-residual [mm]; Hits");
	combination.add(correctedHist, "z-residuals with time walk correction");
//	combination.normalise();
	combination.setNcolumn(1);
	combination.createCombined();
	TGaxis::SetMaxDigits(4);
}

TH1* plotFittedTimeWalk(TH2* uncorrected) {
	uncorrected->FitSlicesY();
	auto means=dynamic_cast<TH1*>( gDirectory->Get("zResidualByToT_1") );
	if(!means) { cerr<<"failed to retrieve result from fitslicesy()\n"; return nullptr; };

	//rebin means
	if (false) {
		auto hist=means;
		std::vector<double> lowerEdges, contents, errors;
		lowerEdges.push_back( hist->GetXaxis()->GetBinLowEdge(1));
		for(int bin=1; bin<=hist->GetNbinsX();) {			
			double lowEdge=hist->GetXaxis()->GetBinLowEdge(bin);
			int nBinsPerBin=  lowEdge<1.5 ? 1 : lowEdge<2 ? 2 : 4;

			double content=0, error2=0;
			for(int i=0; i<nBinsPerBin; i++) {
				content+=hist->GetBinContent(bin);
				error2+=hist->GetBinError(bin)*hist->GetBinError(bin);
				bin++;
			}
			contents.push_back(content/nBinsPerBin);
			errors.push_back(sqrt(error2)/nBinsPerBin);
			lowerEdges.push_back(hist->GetXaxis()->GetBinLowEdge(bin));

		}
		auto rebinned=new TH1D("rebinned", "rebinned", lowerEdges.size()-1, lowerEdges.data() );
		for(int i=0; i<contents.size(); i++) {
			rebinned->SetBinContent(i+1,contents[i]);
			rebinned->SetBinError(i+1,errors[i]);
		}
		means=rebinned;
	}
	

	double minToT=0.10;
	auto simple=new TF1("2 parameters", "[c_{1}]/(x+[t_{0}])+[offset]", minToT,2.5);
	simple->SetLineWidth(1);
	
	simple->SetParameters(0.366,0.066,0);
	means->Fit(simple, "", "", minToT,2.5);
	double offset=simple->GetParameter("offset");
	means=shiftY(means, -simple->GetParameter("offset") );

	new TCanvas("fittedTW", "fitted TW", 800,600);
	means->Draw();	
	simple->SetParameters(0.366,0.066,0);
	means->Fit(simple, "", "", minToT,2.5);


	means->GetYaxis()->SetTitle("Mean z-residual [mm]");
	means->GetYaxis()->SetTitleOffset(1.3);
	means->GetYaxis()->SetTitleSize(0.05);
	means->GetYaxis()->SetLabelSize(0.05);
	means->GetYaxis()->SetRangeUser(0,3);
	means->GetXaxis()->SetTitleOffset(1.3);
	means->GetXaxis()->SetTitleSize(0.05);
	means->GetXaxis()->SetLabelSize(0.05);
	gPad->SetMargin(0.15,0.1,0.15,0.1);

	simple->SetLineWidth(1);
	simple->SetNpx(1000);
	gPad->SetTicks(1,0);
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	gStyle->SetOptTitle(0);

  TGaxis *rightAxis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),gPad->GetUxmax(), gPad->GetUymax(),0,3./0.0546 /*residual divided by drift speed*/,510, "+L");
  rightAxis->SetTitle("[ns]");
  rightAxis->SetTitleOffset(1.0);
  rightAxis->SetTitleSize(0.05);
  rightAxis->SetLabelSize(0.05);
  rightAxis->SetTextFont(42);
  rightAxis->SetLabelFont(42);
  rightAxis->Draw();

	auto stats=new StatsWrapper(0.5,0.6,0.88,0.88);
	stats->add("c_{1}", simple->GetParameter("c_{1}"), 3, "mm #mus");
	stats->add("t_{0}", simple->GetParameter("t_{0}"), 4, "#mus");
//	stats->add("z_{offset}", offset, 3, "mm");
	//stats->addChiSquare(*simple);
	stats->draw();
	return means;
}
TH1* plotFittedTimeWalk( std::string filename="combinedFit_fiducialAll.root", std::string object="quad/zResidualByToT") {
	auto uncorrected=getObjectFromFile<TH2>(object, filename);
	return plotFittedTimeWalk(uncorrected);
}
//draw four timewalk on 4 pads in one canvas
void combineFittedTimeWalkForChips(std::string filename="combinedFit.root", std::string object="zResidualByToT") {
	auto canv=new TCanvas("canvas", "combined canvas", 1000,1000);
	canv->Divide(2,2);
	HistogramCombiner comb{"combinedTW"};
	for(int i=0;i<4;i++) {
		canv->cd(i+1);
		auto h=plotFittedTimeWalk(filename, chipDirectories[i]+"/"+object);
		//comb.add(h, "chip "+to_string(i+1) );
	}
	//comb.createCombined();
}


void combineTimeWalkResiduals(std::string filename="combinedFit.root", std::string object="zResidualByToTCorrected") {
	//auto canv=new TCanvas("canvas", "combined canvas", 1000,1000);
	//canv->Divide(2,2);
	HistogramCombiner comb{"combinedTWResiduals"};
	for(int i=0;i<4;i++) {
		//canv->cd(i+1);
		auto corrected=getObjectFromFile<TH2>(chipDirectories[i]+"/"+object, filename);
		corrected->FitSlicesY();
		auto means=dynamic_cast<TH1*>( gDirectory->Get("zResidualByToTCorrected_1") );
		comb.add(means, "chip "+to_string(i+1) );
	}
	comb.createCombined();
}


//fitSlicesY with specification of range
TH1* fitDiffusionSlices(TH2* h2, std::string x="z") {
	TF1* gaus=new TF1("gaus","gaus(0)", -2,2);
	TF1* gausBG=new TF1("gaus","gaus(0)+[3]", -2,2);
	TF1* gausRange=new TF1("gausRange","gaus(0)", -0.5,0.5);
	gausRange->SetParameters(4E4,0.05,0.22);
//	TF1* exGaus=new TF1("exGaus", "[c]*[l]/2*exp([l]/2*(2*[m]+[l]*[s]*[s]-2*x))*TMath::Erfc( ([m]+[l]*[s]*[s]-x)/sqrt(2)/[s] )", -2, 2); //Exponentially modified gaussian distribution (see wiki)
//	exGaus->SetParameters(2E4,3.1,-0.25,0.2); // Constant, Lambda, Mean, Sigma
	//exGaus->FixParameter(1,3.1);
	gausBG->SetParameters(1E4,0,0.3,0); // Constant, Mean, Sigma, bg
	gausBG->SetParLimits(2, 1E-3,10);
	gaus->SetParameters(1E4,0,0.38); // Constant, Mean, Sigma, bg

	h2->FitSlicesY(x=="z" ? gausRange : gaus, 0/*firstbin*/, -1/*lastbin*/, 15/*min number of entries*/, "QNRM");
//	h2->FitSlicesY(gaus, 0, -1, 50, "QNR");

	//view background contributions
//	gDirectory->Get( (h2->GetName()+std::string("_3")).c_str() )->DrawClone();
//	new TCanvas();

	return dynamic_cast<TH1*>(gDirectory->Get( (h2->GetName()+std::string("_2")).c_str() ));
}


TF1* fitDiffusion( TH2* h2 , std::string x="x", double z0=-1, std::string canvname="canv") {
	double zmax=9;
	std::string LorT=x=="x" ? "T" : x=="z" ? "L" : x;

	TF1* drift=new TF1("drift", ("sqrt( pow([#sigma_{"+x+"0}],2) + pow([D_{"+x+"}],2)/10*(x-[z0]) )").c_str(), z0, zmax);
	drift->SetLineWidth(1);
	new TCanvas((canvname+"_"+x).c_str(), (canvname+"_"+x).c_str(), 800,600);

	fitDiffusionSlices(h2,x);

	TH1* h2_2=nullptr;
	if(false && x=="z") {
		TH1 *h22 =dynamic_cast<TH1*>(gDirectory->Get( (h2->GetName()+std::string("_1")).c_str() ));
		TH1 *h23 =dynamic_cast<TH1*>(gDirectory->Get( (h2->GetName()+std::string("_3")).c_str() ));
		h2_2=h23;//makeOperated(h22, h23, [](double l, double s) { return fabs(l) < 1E-30 ? 0 : sqrt( 1./(l*l)+s*s ); });
	} else {
		h2_2 =dynamic_cast<TH1*>(gDirectory->Get( (h2->GetName()+std::string("_2")).c_str() ));
	}
	if(!h2_2) { auto msg="could not get results from fit\n"; std::cerr<<msg; throw msg; }

	h2_2->GetXaxis()->SetTitle("z-position [mm]");
	h2_2->GetYaxis()->SetTitle( ("#sigma_{"+x+"} from fit to track-residual [mm]").c_str() );
	increaseAxisSize(h2_2, 0.05);	
	h2_2->GetYaxis()->SetRangeUser(0, x=="x"?0.4:0.5 );
	h2_2->GetXaxis()->SetRangeUser(z0,zmax);
	
	//guess parameters
	if(x=="z") drift->FixParameter(2,z0);
	else drift->SetParameter(2,z0); //z0
	drift->SetParameter(1,0.3); //D
	if(x=="x") 	drift->FixParameter(0, 0.0158771);
	else drift->SetParameter(0, 0.15);//sigma0


	//add error because of guard
	if(x=="x")
	{
		for(int i=12; i<=15; i++) {
			h2_2->SetBinError(i, 0.01);
		}
	}

	//add error to fit
//	h2_2=addErrorToHist(h2_2, 1E-3); //set all error bins equal
	h2_2->Fit(drift, "", "", z0, zmax);

	gStyle->SetOptTitle(0);	
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);

	h2_2->Draw();
	gPad->SetTicks(1,1);
	gPad->SetMargin(0.15,0.1,0.15,0.1);

	//nicer stat pane
	gStyle->SetOptFit(false);
	auto stats=new StatsWrapper();
	stats->add("D_{"+LorT+"}",  drift->GetParameter(1)*1E3, 0, "#mum/#sqrt{cm}" );
	if(x!="x") stats->add("#sigma_{"+x+"0}" ,  drift->GetParameter(0)*1E3, 0, "#mum" );
	stats->add("z0" ,  drift->GetParameter(2), 2, "mm" );
	//stats->addChiSquare(*drift);
	stats->draw();
	
	h2_2->GetXaxis()->SetRangeUser(drift->GetParameter(2),zmax);

	//plot residuals
	const bool plotResiduals=false;
	if(plotResiduals) {
		auto residualHistogram=getResidualHistogram(h2_2, drift);
		new TCanvas(("residuals"+x).c_str(), "Residuals");
		residualHistogram->Draw();
	}
	
	return drift;

}


void plotDiffusionFromHist(std::string filename="combinedFit.root", std::string histogramName="locExp/xResidualByz", std::string dir="x") {
	auto hists=chipDirectories;
	hists.push_back("quad");
	for(std::string chip : hists )
//	std::string chip="quad";
	{
		TH2* h2=getObjectFromFile<TH2>( chip+"/"+histogramName, filename);
		fitDiffusion(h2, dir, -0.8, chip);
	}
}


void plotDiffusionUsingCut(
		std::string filename="combinedFit.root",
		std::string histogramName="locExp/xResidualByz",
		std::string dir="x", std::string alignFile="align.dat") {
//	auto hists=chipDirectories;
//	hists.push_back("quad");
//	for(std::string chip : hists )
	std::string chip="quad";
	Alignment align(alignFile);
	{
		TH2* h2=getObjectFromFile<TH2>( chip+"/"+histogramName, filename);
		h2=removeBinsByPosition(h2, [&](double x, double y) {
			return fabs(y) < 2.5*align.hitErrors.hitError(x).X();
		});
		new TCanvas();
		h2->Draw("colz");
		fitDiffusion(h2, dir, -0.8, chip);
	}
}

void plotDiffusionCombined(std::string filename="combinedFit.root") {
	TH2D* histogram=nullptr;
	for(std::string x : {"x", "z"}) {
		for(std::string chip : chipDirectories ) {
			TH2D* h=getObjectFromFile<TH2D>( chip+"/locExp"+x+"ResidualByz", filename);
			if(histogram) histogram->Add(h);
			else histogram=h;
		}
		fitDiffusion(histogram, x, 0);
	}
}

void plotZResidualsByToT(std::string filename="combinedFit.root", std::string histname="zResidualByToTCorrected") {
	HistogramCombiner combination("ToTCombination");
	TF1* gaus=new TF1("gaus","gaus(0)", -2,2);
	for(std::string chip : chipDirectories ) {
		TH2* h2=getObjectFromFile<TH2>( chip+"/"+histname, filename);

		h2->FitSlicesY(gaus, 0/*firstbin*/, -1/*lastbin*/, 5/*min number of entries*/, "QNR");

		auto h1=dynamic_cast<TH1*>(gDirectory->Get( (h2->GetName()+std::string("_1")).c_str() ));
		combination.add(h1, chip);
	}
	combination.createCombined();
}


TH1* plotSpot(std::string filename="fitted.root") {
	auto tree=getObjectFromFile<TTree>("fitResults",filename);
	new TCanvas("spot", "spot", 800,800);
	auto histForCounting=getHistFromTree(tree, "hitAverage.x", "hitAverage.x>0", "count", "");
	int nEntries=histForCounting->GetEntries();
	auto hist=getHistFromTree(tree, "residual.y/0.055:residual.x/0.055", "hitAverage.x>0 && flag>0", "spot(80,-40,40,80,-40,40)", "colz");
	hist->GetXaxis()->SetTitle("Columns");
	hist->GetYaxis()->SetTitle("Rows");
	hist->Scale(1./nEntries);
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	return hist;
}


void plotSlicedDiffusionWithFit( std::string filename="combinedFit.root", double z0=-0.67, std::string object="quad/locExp/zResidualByzByToT" ) {
	auto hist3=getObjectFromFile<TH3D>(object, filename);
	auto zaxis=hist3->GetZaxis();


	//HistogramCombiner slices("slices");
	std::vector<std::string> slicename = {"0.15 #mus < ToT < 0.50 #mus", "ToT > 0.50 #mus"};
	std::vector<std::pair<double,double> > binRanges = { {7,21}, {22, zaxis->GetNbins()+1} }; //0.025 per bin!
	std::vector<TH1*> histograms= { nullptr, nullptr };
	auto stats=new StatsWrapper();
	TF1* drift=new TF1("drift", "sqrt( pow([#sigma_{z0}],2) + pow([D],2)/10*(x-[z0]) )", -1, 10);
	drift->SetLineWidth(1);
	TLegend* legend = new TLegend();
	for(int i : {0,1} ) {
		zaxis->SetRange(binRanges[i].first, binRanges[i].second);
		auto proj=(TH2*)hist3->Project3D("yx");

		auto hsigma=fitDiffusionSlices(proj, "z");

		//add error because of guard?
		for(int i=12; i<=15; i++) {
			hsigma->SetBinError(i, 0.01);
		}

		auto name="slice "+std::to_string(i);
		histograms.at(i)=(TH1*)hsigma->Clone(name.c_str());

		histograms.at(i)->GetXaxis()->SetTitle("z-position [mm]");
		histograms.at(i)->GetYaxis()->SetTitle("#sigma_{z} from fit to track-residual [mm]" );
		//increaseAxisSize(histograms.at(i), 0.05);
		//histograms.at(i)->GetYaxis()->SetRangeUser(0,0.45);
		//histograms.at(i)->GetXaxis()->SetRangeUser(z0,23);

		legend->AddEntry(histograms.at(i), (slicename.at(i)+" ("+to_string_precision(proj->GetEntries()/hist3->GetEntries()*100,0)+"%)" ).c_str() );

		cout<<proj->GetEntries()<<" entries\n";

		//slices.add(clone,slicename.at(i)+";z-position [mm];#sigma_{z} from fit to track-residual [mm]");//+" ("+std::to_string(int(proj->GetEntries()/1000000+0.49))+"M hits)");
		if(i==1) {
			auto driftParams=fitDiffusion(proj, "z", z0);
			//cin.get();
			stats->add( "D_{L}",  driftParams->GetParameter(1)*1E3, 0, "#mum/#sqrt{cm}" );
			stats->add( "#sigma_{z0}" ,  driftParams->GetParameter(0)*1E3, 0, "#mum" );
			stats->add( "z0 (fixed)" ,  driftParams->GetParameter(2), 2, "mm" );

			for(int j=0; j<4; j++) {
				drift->SetParameter(j, driftParams->GetParameter(j));
			}
		}
	}

	changeLegendStyle(legend, 1, 0.055);
	legend->SetX1NDC(0.18);

	//new TCanvas();
	//drift->Draw();
	histograms[0]->SetLineColor(kGreen+2);
	histograms[0]->Draw("same");
	histograms[1]->Draw("same");
	drift->Draw("same");
	legend->Draw();

	//slices.setStyle(10);slices.titleSize=0.05;
//	slices.setYRange({0,0.6});
//	slices.createCombined();

	stats->draw();
	//gPad->SetMargin(0.15,0.1,0.15,0.1);

}

void combineFitDiffusion(
		std::string filename="combinedFit.root",
		std::vector<std::string> slicename = {"diffusion"},
		std::vector<std::string> objects = {"quad/locExp/xResidualByz"},
		double z0=-0.73 ) {

	std::vector<TH1*> histograms{objects.size(), nullptr};
	auto stats=new StatsWrapper();
	TF1* drift=new TF1("drift", "sqrt( pow([#sigma_{z0}],2) + pow([D],2)/10*(x-[z0]) )", -1, 10);
	drift->SetLineWidth(1);
	TLegend* legend = new TLegend();

	for(	int i=0; i<objects.size(); i++ ) {
		auto proj=getObjectFromFile<TH2D>(objects.at(i), filename);


		auto hsigma=fitDiffusionSlices(proj, "x");

		auto name="slice "+std::to_string(i);
		histograms.at(i)=(TH1*)hsigma->Clone(name.c_str());

		histograms.at(i)->GetXaxis()->SetTitle("z-position [mm]");
		histograms.at(i)->GetYaxis()->SetTitle("#sigma_{z} from fit to track-residual [mm]" );
		//increaseAxisSize(histograms.at(i), 0.05);
		//histograms.at(i)->GetYaxis()->SetRangeUser(0,0.45);
		//histograms.at(i)->GetXaxis()->SetRangeUser(z0,23);

		legend->AddEntry(histograms.at(i), (slicename.at(i)).c_str() );



	}

	if(objects.size()==2) {
		auto subtracted=makeOperated(histograms[0], histograms[1], [](double a,double b){return 2*b-a;});
		subtracted->Draw();
		histograms.push_back(subtracted);
		legend->AddEntry(subtracted, "Subtracted");

		drift->SetParameter(2,z0);
		drift->FixParameter(0, 0.0158771);
		subtracted->Fit(drift, "", "", z0, 9);

		gStyle->SetOptTitle(0);
		gStyle->SetOptStat(0);
		gStyle->SetOptFit(1);

		subtracted->Draw();
		gPad->SetTicks(1,1);
		gPad->SetMargin(0.15,0.1,0.15,0.1);

		//nicer stat pane
		gStyle->SetOptFit(false);
		//cin.get();
		stats->add( "D_{L}",  drift->GetParameter(1)*1E3, 0, "#mum/#sqrt{cm}" );
		stats->add( "#sigma_{z0}" ,  drift->GetParameter(0)*1E3, 0, "#mum" );
		stats->add( "z0" ,  drift->GetParameter(2), 2, "mm" );
	}

	changeLegendStyle(legend, 1, 0.055);
	legend->SetX1NDC(0.18);

	//new TCanvas();
	//drift->Draw();
//	histograms[0]->SetLineColor(kGreen+2);
	stats->draw();
	for(int i=0; i<histograms.size(); i++) {
		histograms[i]->Draw("same");
	}
	drift->Draw("same");
	legend->Draw();

}


void plotDriftVelocity( std::string filename="fitted.root", std::string alignFile="align.dat") {
	Alignment align(alignFile);
	auto file=openFile(filename);
	align.updateDriftSpeed(*file);
}

/*
void plotDriftVelocity( std::string filename="fitted.root", std::string alignFile="../align.dat") {
	auto tree=getObjectFromFile<TTree>("fitResults",filename);
	Alignment	alignment(alignFile);
	std::cout<<std::to_string(alignment.driftSpeed.value)<<"\n";
	auto vDrift=getHistFromTree(tree, "hitAverage.z/"+std::to_string(alignment.driftSpeed.value)+":laser.z", "fabs(hitAverage.x)+fabs(hitAverage.y)>0 && nHitsPassed>1", "vDrift(26,-0.5,25.5)", "prof");
	vDrift->SetLineWidth(0);
	vDrift->SetMarkerStyle(7);

	auto driftRelation=new TF1("driftRelation", "x/[vdrift]+[intercept]");
	driftRelation->SetParameters(1,1);
	vDrift->Fit(driftRelation);		

	vDrift->GetXaxis()->SetTitle("Laser position [mm]");
	vDrift->GetYaxis()->SetTitle("Average time of arrival [ns]");
	increaseAxisSize(vDrift,0.05);

	gPad->SetMargin(0.15,0.1,0.15,0.1);
	driftRelation->SetLineWidth(1);
	gPad->SetTicks(1,1);
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	gStyle->SetOptTitle(0);

	auto stats=new StatsWrapper();
	stats->add("v_{drift}", driftRelation->GetParameter("vdrift")*1E3, 1, "#mum/ns");
	stats->draw();
}*/

void testTransformBackFunction() {
	Alignment align{"../align.dat"};
	std::string fileName="fitted.root";
	TFile file(fileName.c_str(), "READ");
	TTreeReader reader("fitResults", &file);
	TTreeReaderArray<double> hitpositionx{reader, "hits.position.x"}, hitpositiony{reader, "hits.position.y"};
	
	while(std::cin.get()!='q' and reader.Next() ) {
		for(int i=0; i<hitpositionx.GetSize(); i++) {
				cout<<hitpositionx[i]<<" "<<hitpositiony[i]<<" --> ";
				cout<<align.transformBack( {hitpositionx[i], hitpositiony[i], 0} )<<"\n";
		}
	}

}

void plotCorrelations( std::string combinedFileName="combinedFit.root") {
	auto fitResults=getObjectFromFile<TTree>("fitResults", combinedFileName);
	new TCanvas();
	fitResults->Draw("XZ.intercept+XZ.slope*172:x", "fabs(YZ.intercept+YZ.slope*172-z)<3");
	new TCanvas();
	fitResults->Draw("YZ.intercept+YZ.slope*172:z", "fabs(XZ.intercept+XZ.slope*172-x)<3");

}


TProfile2D* getCombinedDeformations(std::string histName="quad/locExp/xResidualByPosition",
		std::vector<std::string> fileNames={"./run668/combinedFit.root","./run672/combinedFit.root","./run676/combinedFit.root"}) {
	auto prof=getObjectFromFile<TProfile2D>(histName, fileNames.front());
	for(int i=1; i<fileNames.size(); i++) {
		auto p=getObjectFromFile<TProfile2D>(histName, fileNames[i]);
		prof->Add(p);
	}

	return prof;
}

TH2D* getCombinedHitMap(std::string histName="quad/locExp/positionHitMap",
		std::vector<std::string> fileNames={"./run668/combinedFit.root","./run672/combinedFit.root","./run676/combinedFit.root"}) {
	auto prof=getObjectFromFile<TH2D>(histName, fileNames.front());
	for(int i=1; i<fileNames.size(); i++) {
		auto p=getObjectFromFile<TH2D>(histName, fileNames[i]);
		prof->Add(p);
	}
	return prof;
}

TCanvas* plotCombinedHitMap( TH2D* hist ) {
	auto canv=new TCanvas( (hist->GetName()+std::string("deformcanv")).c_str(), "deformations", 900,1000)	;

	hist->SetTitle("");
	hist->GetYaxis()->SetTitleOffset(1.3);
	hist->GetXaxis()->SetTitleOffset(1.1);
	hist->GetZaxis()->SetTitleOffset(1.6);
	hist->Draw("colz0");

	//gStyle->SetPalette(kRainBow);
	gPad->SetMargin(0.1,0.2,0.1,0.05);
	gPad->SetTicks(1,1);

	drawChipEdgesLocal("./run668/align.dat");
	return canv;
}

std::vector<std::vector<double>> fitDeformationCorrections(
		std::vector<std::string> fileNames={"./run668/combinedFit.root","./run672/combinedFit.root","./run676/combinedFit.root"},
		std::string histName="locExp/xResidualByPosition"
) {
	auto correction=new TF1("correction",
			"breitwigner(0)+breitwigner(3) "
			"+"
			"breitwigner(6)+breitwigner(9)+[offset]");
	std::vector<std::vector<double> > estimatedValues={
			{7,1,2.6, /*breit wigner 1*/ -3,2,3, /*breit wigner 1*/3, 14, 3, /*breit wigner 4*/ -8, 15, 2,/*breit wigner 4*/ 0 /*offset*/ },
			{7,1,2.6, /*breit wigner 1*/ -3,2,3, /*breit wigner 1*/3, 14, 3, /*breit wigner 4*/ -8, 15, 2,/*breit wigner 4*/ 0 /*offset*/ },
			{7,14,2.6, /*breit wigner 1*/ -3,15,3, /*breit wigner 1*/3, 27, 3, /*breit wigner 4*/ -8, 28, 2,/*breit wigner 4*/ 0 /*offset*/ },
			{7,14,2.6, /*breit wigner 1*/ -3,15,3, /*breit wigner 1*/3, 27, 3, /*breit wigner 4*/ -8, 28, 2,/*breit wigner 4*/ 0 /*offset*/ } };

	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);

	auto canv=new TCanvas();
	canv->Divide(2,2);

	std::cout<<"{";
	for(int i=0; i<4; i++) {
		canv->cd(i+1);
		gPad->SetTicks(1,1);

//		auto xResidualByPosition=getObjectFromFile<TProfile2D>(chipDirectories[i]+"/xResidualByPosition_locExp", fileName );
		auto xResidualByPosition=getCombinedDeformations(chipDirectories[i]+"/"+histName,fileNames);

		auto xResidualByx=xResidualByPosition->ProfileX();

		for(int j=0; j<correction->GetNpar(); j++) correction->SetParameter( j, estimatedValues[i][j] );

		xResidualByx->Fit(correction, "QS", "");//, chipRange[i].min, chipRange[i].max);
		xResidualByx->SetTitle(";x-position [mm];x-residual [mm]");

		std::cout<<"{";
		for(int j=0; j<correction->GetNpar(); j++) {
			if(j==correction->GetNpar()-1) std::cout<<correction->GetParameter(j)<<( i!=3 ? "},\n" : "}}\n");
			else std::cout<<correction->GetParameter(j)<<", ";
			estimatedValues[i][j]=correction->GetParameter(j);
		}
	}

	return estimatedValues;
}



void fitDeformationCorrectionsPerSliceZ(
		std::vector<std::string> fileNames={"./run668/combinedFit.root","./run672/combinedFit.root","./run676/combinedFit.root"}) {

	auto correction=new TF1("correction",
			"breitwigner(0)+breitwigner(3) "
			"+"
			"breitwigner(6)+breitwigner(9)+[offset]");
	std::vector<std::vector<double> > estimatedValues={
			{7,1,2.6, /*breit wigner 1*/ -3,2,3, /*breit wigner 1*/3, 14, 3, /*breit wigner 4*/ -8, 15, 2,/*breit wigner 4*/ 0 /*offset*/ },
			{7,1,2.6, /*breit wigner 1*/ -3,2,3, /*breit wigner 1*/3, 14, 3, /*breit wigner 4*/ -8, 15, 2,/*breit wigner 4*/ 0 /*offset*/ },
			{7,14,2.6, /*breit wigner 1*/ -3,15,3, /*breit wigner 1*/3, 27, 3, /*breit wigner 4*/ -8, 28, 2,/*breit wigner 4*/ 0 /*offset*/ },
			{7,14,2.6, /*breit wigner 1*/ -3,15,3, /*breit wigner 1*/3, 27, 3, /*breit wigner 4*/ -8, 28, 2,/*breit wigner 4*/ 0 /*offset*/ } };

	const std::string parNames[13] = { "scale", "mean", "sigma", "scale", "mean", "sigma", "scale", "mean", "sigma", "scale", "mean", "sigma", "offset" };


	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);

	auto xByXZCanvas=new TCanvas("xByXZ", "xByXZ", 1000,1000);
	xByXZCanvas->Divide(2,2);
	auto fitParametersCanvas=new TCanvas("fitSlices", "fitSlices", 1000,1000);
	fitParametersCanvas->Divide(2,2);

	std::cout<<"{";
	for(int iChip=0; iChip<4; iChip++) {
		xByXZCanvas->cd(iChip+1);
		gPad->SetTicks(1,1);
		auto xResidualByXZ=getCombinedDeformations(chipDirectories[iChip]+"/locExp/xResidualByXZ",fileNames);
		xResidualByXZ->Rebin2D(4,4);
		xResidualByXZ->SetAxisRange(-2,10,"Y");
		xResidualByXZ->Draw("colz");

		std::vector<TH1*> fittedParams;
		for(int j=0; j<13; j++) {
			auto h=newProjectionHistogram(xResidualByXZ, "_c"+std::to_string(iChip)+"f"+std::to_string(j), "x");
			h->SetTitle(("Chip "+to_string(iChip)+";z-position [mm]; parameter '"+parNames[j]+"'" ).c_str() );
			fittedParams.push_back( h );
		}

		fitParametersCanvas->cd(iChip+1);
		int iSlice=0;
		forEachSlice(xResidualByXZ, "x", [&](TH1* proj) {
			iSlice++;
			if(proj->GetEntries()<1000) return;
			for(int j=0; j<correction->GetNpar(); j++) correction->SetParameter( j, estimatedValues[iChip][j] );
			proj->Fit(correction, "QS", "");//, chipRange[i].min, chipRange[i].max);

			for(int j=0; j<correction->GetNpar()-1; j++) {
				fittedParams[j]->SetBinContent(iSlice, correction->GetParameter(j));
				fittedParams[j]->SetBinError( iSlice, correction->GetParError(j));
			}

//			proj->DrawCopy();
//			gPad->Update();
//			std::cin.get();
			return;
		});

		HistogramCombiner comb("comb", {fittedParams[0], fittedParams[3], fittedParams[6], fittedParams[9]});
		comb.setStyle(7);
		comb.makeLegend=false;
		comb.createCombinedOnCanvas(gPad);
	}
}

std::vector<std::vector<double>> fitDeformationCorrectionsPerSliceY(
		std::vector<std::string> fileNames={"./run668/combinedFit.root","./run672/combinedFit.root","./run676/combinedFit.root"}) {

	auto correction=new TF1("correction",
			"breitwigner(0)+breitwigner(3) "
			"+"
			"breitwigner(6)+breitwigner(9)+[offset]");
	std::vector<std::vector<double> > estimatedValues=fitDeformationCorrections(fileNames);

	const std::string parNames[13] = { "scale", "mean", "sigma", "scale", "mean", "sigma", "scale", "mean", "sigma", "scale", "mean", "sigma", "offset" };


	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);

	auto xByXZCanvas=new TCanvas("xByXY", "xByXY", 1000,1000);
	xByXZCanvas->Divide(2,2);
	auto fitParametersCanvas=new TCanvas("fitSlices", "fitSlices", 1000,1000);
	fitParametersCanvas->Divide(2,2);

	for(int iChip=0; iChip<4; iChip++) {
//		std::cout<<"\n\n";
		xByXZCanvas->cd(iChip+1);
		gPad->SetTicks(1,1);
		auto xResidualByXY=getCombinedDeformations(chipDirectories[iChip]+"/locExp/xResidualByPosition",fileNames);
		xResidualByXY->Rebin2D(4,16);
		//xResidualByXZ->SetAxisRange(-2,10,"Y");
		xResidualByXY->Draw("colz0");

		std::vector<TH1*> fittedParams;
		for(int j=0; j<13; j++) {
			auto h=newProjectionHistogram(xResidualByXY, "_c"+std::to_string(iChip)+"f"+std::to_string(j), "x");
			h->SetTitle(("Chip "+to_string(iChip)+";y-position [mm]; parameter '"+parNames[j]+"'" ).c_str() );
			fittedParams.push_back( h );
		}

		const int selectedParameter=0;
		fitParametersCanvas->cd(iChip+1);
		int iSlice=0;
		forEachSlice(xResidualByXY, "x", [&](TH1* proj) {
			iSlice++;
			if(proj->GetEntries()<100) return;
			for(int j=0; j<correction->GetNpar(); j++) {
				if(j%3==selectedParameter)
					correction->SetParameter( j, estimatedValues[iChip][j] );
				else
					correction->FixParameter( j, estimatedValues[iChip][j] );
			}
			proj->Fit(correction, "QS", "");//, chipRange[i].min, chipRange[i].max);

			for(int j=0; j<correction->GetNpar()-1; j++) {
				auto parameterValue=correction->GetParameter(j);
				fittedParams[j]->SetBinContent(iSlice, parameterValue );
				fittedParams[j]->SetBinError( iSlice, correction->GetParError(j));
			}

//			proj->DrawCopy();
//			gPad->Update();
//			std::cin.get();
			return;
		});

		estimatedValues[iChip].resize(33);
		estimatedValues[iChip][32]=estimatedValues[iChip][12];
		for(int iPar=0; iPar<4; iPar++) {
			auto fitResult = fittedParams[iPar*3+selectedParameter]->Fit("pol4", "QSM");
			auto aCoefficient=fitResult->GetParams()[0];
			estimatedValues[iChip][iPar*3+selectedParameter]=aCoefficient;
			for(int iFit=1; iFit<5; iFit++) {
				estimatedValues[iChip][12+iPar+4*(iFit-1)]=fitResult->GetParams()[iFit]/aCoefficient;
//				std::cout<<iFit<<": "<<fitResult->GetParams()[iFit]/aCoefficient<<"\n";
			}
		}

		HistogramCombiner comb("comb", {fittedParams[0+selectedParameter], fittedParams[3+selectedParameter], fittedParams[6+selectedParameter], fittedParams[9+selectedParameter]});
		comb.setStyle(7);
		comb.makeLegend=false;
		comb.createCombinedOnCanvas(gPad);
	}

	return estimatedValues;
}


std::vector<std::vector<double>> fitDeformationCorrections2D(
		std::vector<std::string> fileNames={"./run668/combinedFit.root","./run672/combinedFit.root","./run676/combinedFit.root"},
		std::string histName="locExp/xResidualByPosition",
		std::string alignFile="run668/align.dat"
) {
	auto correction=new TF2("correction",
			"(1+[b0]*y+[c0]*y*y+[d0]*y*y*y+[e0]*y*y*y*y)*breitwigner(0)+(1+[b1]*y+[c1]*y*y+[d1]*y*y*y+[e1]*y*y*y*y)*breitwigner(3) "
			"+"
			"(1+[b2]*y+[c2]*y*y+[d2]*y*y*y+[e2]*y*y*y*y)*breitwigner(6)+(1+[b3]*y+[c3]*y*y+[d3]*y*y*y+[e3]*y*y*y*y)*breitwigner(9)+[offset]");
	correction->SetLineWidth(1);
	correction->SetNpx(1000);
	correction->SetNpy(1000);

	std::vector<std::vector<double> > estimatedValues=fitDeformationCorrectionsPerSliceY(fileNames);
//	{
//			{7,1,2.6, /*breit wigner 1*/ -3,2,3, /*breit wigner 1*/3, 14, 3, /*breit wigner 4*/ -8, 15, 2,/*breit wigner 4*/   1,1,1,1, 1,1,1,1, 1E-2,1E-2,1E-2,1E-2, 1E-3,1E-3,1E-3,1E-3, 0,0,0,0, /*a,b,c,d,e*/ 0 /*offset*/ },
//			{7,1,2.6, /*breit wigner 1*/ -3,2,3, /*breit wigner 1*/3, 14, 3, /*breit wigner 4*/ -8, 15, 2,/*breit wigner 4*/   1,1,1,1, 1,1,1,1, 1E-2,1E-2,1E-2,1E-2, 1E-3,1E-3,1E-3,1E-3, 0,0,0,0, /*a,b*/  0 /*offset*/ },
//			{7,14,2.6, /*breit wigner 1*/ -3,15,3, /*breit wigner 1*/3, 27, 3, /*breit wigner 4*/ -8, 28, 2,/*breit wigner 4*/ 1,1,1,1, 1,1,1,1, 1E-2,1E-2,1E-2,1E-2, 1E-3,1E-3,1E-3,1E-3, 0,0,0,0, /*a,b*/  0 /*offset*/ },
//			{7,14,2.6, /*breit wigner 1*/ -3,15,3, /*breit wigner 1*/3, 27, 3, /*breit wigner 4*/ -8, 28, 2,/*breit wigner 4*/ 1,1,1,1, 1,1,1,1, 1E-2,1E-2,1E-2,1E-2, 1E-3,1E-3,1E-3,1E-3, 0,0,0,0, /*a,b*/  0 /*offset*/ } };

//	auto fitted1DValues = fitDeformationCorrections(fileNames, histName);
//
//	for(int i=0; i<4; i++) {
//		for(int j=0; j<fitted1DValues[i].size()-1; j++) {
//			estimatedValues[i][j]=fitted1DValues[i][j];
//		}
//		estimatedValues[i][32]=fitted1DValues[i][12];
//	}

	Alignment align(alignFile);


	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);

	auto canv=new TCanvas("2dCorrectinoFit","2D fit of deformations",1000,1200);
//	canv->Divide(2,2);
	canv->SetMargin(0.1,0.2,0.1,0.05);

	std::cout<<"{";
	for(int i=0; i<4; i++) {
//		canv->cd(i+1);
		gPad->SetTicks(1,1);

		auto corners=align.chips[i].getChipCorners();
		auto xminmax=std::minmax_element(corners.begin(), corners.end(), [](const TVector3& a, const TVector3& b){ return a.x() < b.x(); });
		auto yminmax=std::minmax_element(corners.begin(), corners.end(), [](const TVector3& a, const TVector3& b){ return a.y() < b.y(); });

//		auto xResidualByPosition=getObjectFromFile<TProfile2D>(chipDirectories[i]+"/xResidualByPosition_locExp", fileName );
		auto xResidualByPosition=getCombinedDeformations(chipDirectories[i]+"/"+histName,fileNames);
		xResidualByPosition->Rebin2D(4,4);
		removeBinsWithFewerEntries(xResidualByPosition, 800);
//		xResidualByPosition->SetAxisRange( xminmax.first->x(),  xminmax.second->x() + (i==1 ? -0.5 : 0), "X");
//		xResidualByPosition->SetAxisRange( yminmax.first->y(),  yminmax.second->y(), "Y");
//		std::cout<<"\n"<<"yrange="<<yminmax.first->y()<<" - "<<yminmax.second->y()<<"\n";
		xResidualByPosition->SetTitle(";x-position [mm]; y-position [mm]; mean x-residual [mm]");
		xResidualByPosition->SetTitleOffset(1.2, "y");
		xResidualByPosition->SetTitleOffset(1.7, "z");
		xResidualByPosition->Draw(i ? "col0same" : "colz0");

		std::cout<<"/*from slices, not fitted:*/";
		for(int j=0; j<correction->GetNpar(); j++) correction->FixParameter( j, estimatedValues[i][j] );

		//const double margin[]={-0.1,-0.1,-0.05,0};
		const double margin[]={0.5,0.5,0.5,0.5};
		correction->SetRange( xminmax.first->x()+margin[i], yminmax.first->y()+0.3*margin[i], xminmax.second->x()-margin[i], yminmax.second->y()-0.3*margin[i] );
		correction->SetContour(10);

		xResidualByPosition->Fit(correction, "QSR", "");//, chipRange[i].min, chipRange[i].max);
		//xResidualByPosition->SetTitle(";x-position [mm];x-residual [mm]");

		std::cout<<"{";
		for(int j=0; j<correction->GetNpar(); j++) {
			if(j==correction->GetNpar()-1) std::cout<<correction->GetParameter(j)<<( i!=3 ? "},\n" : "}}\n");
			else std::cout<<correction->GetParameter(j)<<", ";
			estimatedValues[i][j]=correction->GetParameter(j);
		}
	}

	align.drawChipEdges(false);

	return estimatedValues;
}


void compareDeformations(std::vector<std::string> fileNames={"./run668/combinedFit.root","./run672/combinedFit.root","./run676/combinedFit.root"}) {

	std::string dir("quad");
//	for(const auto& dir :chipDirectories)
	{
	//	auto xResidualByPosition=getObjectFromFile<TProfile2D>("quad/xResidualByPosition_locExp", fileName );
		auto xResidualByPosition=getCombinedDeformations(dir+"/locExp/xResidualByPosition",fileNames );
		auto xResidualByx=xResidualByPosition->ProfileX();

	//	auto xResidualByPosition_cor=getObjectFromFile<TProfile2D>("quad/xResidualByPosition_corrected", fileName );
		auto xResidualByPosition_cor=getCombinedDeformations(dir+"/corrected/xResidualByPosition",fileNames);
		auto xResidualByx_cor=xResidualByPosition_cor->ProfileX();

		HistogramCombiner deform("deform"+dir);
		deform.add(xResidualByx, "Before correction;x-position [mm]; x-residual [mm]");
		deform.add(xResidualByx_cor, "After correction");
		deform.createCombined();
	}
}


//In QUAD Frame
std::vector< std::vector<TVector3> > getChipAreaFromEdge(double dist=2, std::string alignFile="align.dat") {
	std::vector< std::vector<TVector3> > areas;
	Alignment align(alignFile);
	for(int i=0; i<4; i++) {
				auto corners=align.chips[i].getChipCorners();
				TVector3 meanPos=std::accumulate(corners.begin(), corners.end(), TVector3());
				meanPos*=1.0/corners.size();
				areas.push_back({});
				for(const auto& c : corners) {
					TVector3 unit=(c-meanPos).Unit();
					areas.back().emplace_back( c-dist*unit );
				}
	}
	return areas;
}

void drawAreaPolyLine(const std::vector<TVector3>& corners, Color_t color=kBlack) {
	TPolyLine l;
	for(auto& corner : corners) {
		l.SetNextPoint(corner.x(), corner.y());
	}
	l.SetNextPoint(corners[0].x(), corners[0].y());
	l.SetLineColor(color);
	l.DrawClone();
}

//do frequency on an unweighted profile, i.e. all entreis should have same weight
TH1D* getFrequencyInAreaFromProfile(TProfile2D* original, std::vector<std::vector<TVector3>> areas, double min=-0.1, double max=0.1, int nBins=80, double entryWeight=1.0) {
	TH1D* frequencyHist=new TH1D(
		(original->GetName()+std::string("freq")).c_str(),
		(original->GetTitle()+std::string(" frequency;")+original->GetZaxis()->GetTitle()+";" ).c_str(),
		nBins, min, max);
	for(int i=1;i<=original->GetNbinsX();i++) {
		for(int j=1; j<=original->GetNbinsY();j++) {
			int bin=original->GetBin(i,j);
			if( original->GetBinEntries(bin) <= 0 ) continue;
			double binContent= original->GetBinContent(i,j); //possible rounding error here
			double binError=original->GetBinError(i,j);
//			frequencyHistAll->Fill(binContent);
			TVector3 binPoint{original->GetXaxis()->GetBinCenter(i), original->GetYaxis()->GetBinCenter(j),0};
			for(auto a : areas ) if( isInArea(binPoint, a) ) {
				frequencyHist->Fill(binContent);
				break;
			}
//				frequencyPull->Fill(binContent/error);
		}
	}
	return frequencyHist;
}
TH1D* getFrequencyInArea(TH2D* original, std::vector<std::vector<TVector3>> areas, double min=-0.1, double max=0.1, int nBins=80, double entryWeight=1.0) {
	TH1D* frequencyHist=new TH1D(
		(original->GetName()+std::string("freq")).c_str(),
		(original->GetTitle()+std::string(" frequency;")+original->GetZaxis()->GetTitle()+";" ).c_str(),
		nBins, min, max);
	for(int i=1;i<=original->GetNbinsX();i++) {
		for(int j=1; j<=original->GetNbinsY();j++) {
			int bin=original->GetBin(i,j);
			if( fabs(original->GetBinContent(bin)) < 1E-20 ) continue;
			double binContent= original->GetBinContent(i,j); //possible rounding error here
			double binError=original->GetBinError(i,j);
//			frequencyHistAll->Fill(binContent);
			TVector3 binPoint{original->GetXaxis()->GetBinCenter(i), original->GetYaxis()->GetBinCenter(j),0};
			for(auto a : areas ) if( isInArea(binPoint, a) ) {
				frequencyHist->Fill(binContent);
				break;
			}
//				frequencyPull->Fill(binContent/error);
		}
	}
	return frequencyHist;
}

void drawDeformationsWithAreas(TProfile2D* prof) {
	new TCanvas( (prof->GetName()+std::string("deformcanv")).c_str(), "deformations", 900,1000)	;

	prof->SetTitle("");
	prof->GetYaxis()->SetTitleOffset(1.3);
	prof->GetXaxis()->SetTitleOffset(1.1);
	prof->GetZaxis()->SetTitleOffset(1.6);
	prof->Draw("colz0");

	gStyle->SetPalette(kRainBow);
	gPad->SetMargin(0.1,0.2,0.1,0.05);
	gPad->SetTicks(1,1);

	drawChipEdgesLocal("./run668/align.dat");

	auto areas=getChipAreaFromEdge(3, "./run668/align.dat");
	for(const auto& a : areas) drawAreaPolyLine(a);
	//new TCanvas();
	gStyle->SetOptStat(1);
	HistogramCombiner combined(prof->GetName()+std::string("freq"));
	auto freq=getFrequencyHistogramFromProfile(prof);
	combined.add(freq, "All bins (RMS="+std::to_string( int( freq->GetRMS()*1000) )+" #mum)");
	auto freqArea=getFrequencyInAreaFromProfile(prof, areas);
	combined.add(freqArea, "Bins inside selected area (RMS="+std::to_string( int( freqArea->GetRMS()*1000) )+" #mum)");
	combined.setNcolumn(1);
	combined.createCombined();
}
void drawDeformationsWithAreas(TH2D* h2) {
	new TCanvas( (h2->GetName()+std::string("deformcanv")).c_str(), "deformations", 900,1000)	;

	h2->SetTitle("");
	h2->GetYaxis()->SetTitleOffset(1.3);
	h2->GetXaxis()->SetTitleOffset(1.1);
	h2->GetZaxis()->SetTitleOffset(1.6);
	h2->Draw("colz1");

	gStyle->SetPalette(kRainBow);
	gPad->SetMargin(0.1,0.2,0.1,0.05);
	gPad->SetTicks(1,1);

	drawChipEdgesLocal("./run668/align.dat");

	auto areas=getChipAreaFromEdge(3, "./run668/align.dat");
	for(const auto& a : areas) drawAreaPolyLine(a);
	//new TCanvas();
	gStyle->SetOptStat(1);
	HistogramCombiner combined(h2->GetName()+std::string("freq"));
	auto freq=getFrequencyHistogramFromProfile(h2);
	combined.add(freq, "All bins (RMS="+std::to_string( int( freq->GetRMS()*1000) )+" #mum)");
	auto freqArea=getFrequencyInArea(h2, areas);
	combined.add(freqArea, "Bins inside selected area (RMS="+std::to_string( int( freqArea->GetRMS()*1000) )+" #mum)");
	combined.setNcolumn(1);
	combined.createCombined();
}

TProfile2D* combineDeformations(
		std::string histName="quad/locExp/xResidualByPosition",
		std::vector<std::string> fileNames={"./run668/combinedFit.root","./run672/combinedFit.root","./run676/combinedFit.root"}) {

	auto prof=getCombinedDeformations(histName, fileNames);

	//rebin
	prof->Rebin2D(4,4);
	removeBinsWithFewerEntries(prof, 800);

	drawDeformationsWithAreas(prof);
	return prof;
}

void combineDeformationFrequencies(
		std::vector<std::string> histNames={"quad/locExp/xResidualByPosition",  "quad/corrected/xResidualByPosition" },
		std::vector<std::string> histTitles={"Before correction", "After correction"},
		std::vector<std::string> fileNames={"./run668/combinedFit.root","./run672/combinedFit.root","./run676/combinedFit.root"}) {
	HistogramCombiner combined("combFreq");
	auto histTitle=histTitles.begin();
	for(auto& histName : histNames) {
		auto prof=getCombinedDeformations(histName, fileNames);
		//rebin
		prof->Rebin2D(4,4);
		removeBinsWithFewerEntries(prof, 1000);
		auto hist=getFrequencyHistogramFromProfile(prof);
		combined.add(hist, *histTitle++ + "( RMS = "+to_string( int( hist->GetRMS()*1000) )+" #mum)");

	}
	combined.setNcolumn(1);
	combined.createCombined();
}




TProfile2D* correctDeformation(
		std::vector<std::string> fileNames={"./run668/combinedFit.root","./run672/combinedFit.root","./run676/combinedFit.root"},
		std::string histName="locExp/xResidualByPosition",
		std::string alignFile="run668/align.dat") {
	auto parameters=fitDeformationCorrections(fileNames, histName);
	auto xResidualByPosition=getCombinedDeformations("quad/"+histName,fileNames);

	Alignment align(alignFile);
	auto corrected=makeAppliedByPosition(xResidualByPosition, [&](double xres, double x, double y) {
		auto chipNumber=align.getChipNumber({x,y,0});
		if(chipNumber) xres-=deformationCorrection(chipNumber-1, x, parameters);
		return xres;
	}, "corrected");

	corrected->SetMinimum(-0.1);
	corrected->SetMaximum(0.1);

	//rebin
	corrected->Rebin2D(4,4);
	removeBinsWithFewerEntries(corrected, 800);

	drawDeformationsWithAreas(corrected);
	return corrected;
}


TProfile2D* correctDeformation2D(
		std::vector<std::string> fileNames={"./run668/combinedFit.root","./run672/combinedFit.root","./run676/combinedFit.root"},
		std::string histName="locExp/xResidualByPosition",
		std::string alignFile="run668/align.dat") {
	auto parameters=fitDeformationCorrections2D(fileNames); //fitDeformationCorrectionsPerSliceY(fileNames);
	auto xResidualByPosition=getCombinedDeformations("quad/"+histName,fileNames);

	xResidualByPosition->Rebin2D(4,4);

	Alignment align(alignFile);
	auto corrected=makeAppliedByPosition(xResidualByPosition, [&](double xres, double x, double y) {
		auto chipNumber=align.getChipNumber({x,y,0});
		const double extraOffset=0.0034; //this was found as a result of the slicing fiting procedure, which does not properly fit the offset
		if(chipNumber) xres-=deformationCorrection2D(chipNumber-1, x, y, parameters)+extraOffset;
		return xres;
	}, "corrected");

	corrected->SetMinimum(-0.1);
	corrected->SetMaximum(0.1);

	//rebin
//	corrected->Rebin2D(4,4);
	removeBinsWithFewerEntries(corrected, 800);

	drawDeformationsWithAreas(corrected);
	return corrected;
}


void combineDeformationFrequenciesXandY(
		std::vector<std::string> fileNames={"./run668/combinedFit.root","./run672/combinedFit.root","./run676/combinedFit.root"} ) {

	auto deformx=correctDeformation2D(fileNames);
	auto deformy=combineDeformations("quad/locExp/zResidualByPosition", fileNames);

	const double specialOffset=-0.0044;
	deformy=makeAppliedByPosition(deformy, [specialOffset](double o,double,double) {
		return o-specialOffset;
	});
	drawDeformationsWithAreas(deformy);

	auto areas=getChipAreaFromEdge(3, "./run668/align.dat");

	//new TCanvas();
	gStyle->SetOptStat(1);
	HistogramCombiner combined("combinedFrequenciesXandY");


	auto freqxArea=getFrequencyInArea(deformx, areas);
	std::cout<<"mean="<<freqxArea->GetMean()<<"\n";
	combined.add(freqxArea, "x-residual bins inside selected area (RMS="+std::to_string( int( freqxArea->GetRMS()*1000) )+" #mum)");

	auto freqx=getFrequencyHistogramFromProfile(deformx);
	std::cout<<"mean="<<freqx->GetMean()<<"\n";
	combined.add(freqx, "All x-residual bins (RMS="+std::to_string( int( freqx->GetRMS()*1000) )+" #mum)");

	auto freqyArea=getFrequencyInArea(deformy, areas);
	std::cout<<"mean="<<freqyArea->GetMean()<<"\n";
	combined.add(freqyArea, "z-residual bins inside selected area (RMS="+std::to_string( int( freqyArea->GetRMS()*1000) )+" #mum)");

	auto freqy=getFrequencyHistogramFromProfile(deformy);
	std::cout<<"mean="<<freqy->GetMean()<<"\n";
	combined.add(freqy, "All z-residual bins (RMS="+std::to_string( int( freqy->GetRMS()*1000) )+" #mum)");

	combined.setNcolumn(1);
	combined.setStyle(4);
	combined.createCombined();
}

void compareDeformationSlices(
		std::vector<std::string> fileNames={"./run668/combinedFit.root","./run672/combinedFit.root","./run676/combinedFit.root"},
		std::string name1="locExp/xResidualByXYZ2", std::string name2="locExp/xResidualByXYZ3"
) {
	auto slice2=correctDeformation(fileNames,name1);
	auto slice3=correctDeformation(fileNames,name2);

//	auto difference=(TProfile2D*) slice2->Clone("difference");
//	difference->Add(slice3,-1);
	auto difference=slice2->ProjectionXY("difference");
	difference->Add(slice3->ProjectionXY(),-1);
	removeBinsByPosition(difference, [&](double x, double y) {
		return slice2->GetBinEntries( slice2->FindBin(x,y) ) > 1;
	});
	removeBinsByPosition(difference, [&](double x, double y) {
		return slice3->GetBinEntries( slice3->FindBin(x,y) ) > 1;
	});
	difference->SetMinimum(-0.1);
	difference->SetMaximum(0.1);

	drawDeformationsWithAreas(difference);
}


void plotDeformationPerSlice(std::vector<std::string> fileNames={"./run668/combinedFit.root","./run672/combinedFit.root","./run676/combinedFit.root"}, std::string objectName="quad/locExp/xResidualByXYZ" ) {
	//auto profile3d=getObjectFromFile<TProfile3D>(objectName, fileName);
	//auto zaxis=profile3d->GetZaxis();

	auto canvas=new TCanvas("deformationPerSlice", "deformatationPerSlice", 1200,1800);
	canvas->Divide(3,2);
//	canvas->DivideSquare(zaxis->GetNbins());
	gStyle->SetOptStat(0);

	for(int i=1; i<=5; i++) {
		//zaxis->SetRange(i,i);
		canvas->cd(i);

//		auto prof2d=(TH2*)profile3d->Project3D( ("yx"+std::to_string(i)).c_str() );
//		auto prof2d=getObjectFromFile<TProfile2D>(objectName+std::to_string(i-1), fileName);
		auto prof2d = getCombinedDeformations(objectName+std::to_string(i-1), fileNames);
		prof2d->Rebin2D(8,8);

		removeBinsWithFewerEntries(prof2d, 500);

		prof2d->SetMaximum(0.1);
		prof2d->SetMinimum(-0.1);
		prof2d->Draw("colz0");
	}
}

void plotDeformationsXZ( std::string fileName, std::string objectName="quad/corrected/xResidualByXZ") {
	auto prof2d=getObjectFromFile<TProfile2D>(objectName, fileName);
	prof2d->Rebin2D(6,2);
	removeBinsWithFewerEntries(prof2d, 500);
	gStyle->SetPalette(kRainBow);
	prof2d->SetMinimum(-0.1);
	prof2d->SetMaximum(0.1);
	prof2d->Draw("colz0");
	gPad->SetMargin(0.10,0.15,0.1,0.05);
	gPad->SetTicks(1,1);
	increaseAxisSize(prof2d);
	prof2d->GetZaxis()->SetTitleSize(0.05);
	prof2d->GetZaxis()->SetLabelSize(0.05);
	prof2d->GetYaxis()->SetTitleOffset(1.05);

	HistogramCombiner combined{"combinedxByXZ"};
	auto freq=getFrequencyHistogramFromProfile(prof2d);
	combined.add(freq, "Mean x-residual per bin (RMS="+std::to_string( int( freq->GetRMS()*1000) )+" #mum)");
	combined.setNcolumn(1);
	combined.setStyle(3);
	combined.createCombined();
}

void compareChipEdges(std::vector<std::string> alignFiles) {
	TH2D axis{"axisobj", ";x-position [mm];y-position [mm]",1000,0,29.04,1000,-14.55,25.05 };
	new TCanvas("chipEdges", "chipEdges", 850,1000)	;
	axis.DrawClone();
	const std::vector<Color_t> colors={kBlack,kRed,kBlue,kGreen+1};
	for(int i=0; i<alignFiles.size(); i++) {
		Alignment alignment(alignFiles[i]);
		alignment.drawChipEdges(false, colors[i%colors.size()]);
	}
	Alignment nominal(alignFiles[0]);
	std::cout<<"\n";
	for(int i=1; i<alignFiles.size(); i++){
		std::cout<<alignFiles[i]<<"\n";
		Alignment alignment(alignFiles[i]);
		for(int j=0; j<4; j++) {
			auto corners=alignment.chips[j].getChipCorners();
			auto nominalCorners=nominal.chips[j].getChipCorners();
			std::cout<<"Chip "<<j<<"\n"
			"First pad: "<<1E3*(corners[0].X()-nominalCorners[0].X())<<", "<<1E3*(corners[0].Y()-nominalCorners[0].Y())<<"\n"
			"Second pad: "<<1E3*(corners[1].X()-nominalCorners[1].X())<<", "<<1E3*(corners[1].Y()-nominalCorners[1].Y())<<"\n";
		}	
		std::cout<<"\n";
	}
	gStyle->SetOptStat(0);
}

void drawPartialFitResidual(std::string filename="combinedFit.root") {
	auto tree = getObjectFromFile<TTree>("fitResults",filename);

	for(std::string x : {"x","z"} ) {
		try{
		combineHistogramsFromTree(*tree, { "meanDiffPerChipPerFitLast[0][1]."+x,	"meanDiffPerChipPerFitFirst[0][1]."+x, "meanDiffPerChip[0]."+x},
				"fabs(meanDiffPerChipPerFitLast[0][1]."+x+")>1E-10 && fabs(meanDiffPerChipPerFitFirst[0][1]."+x+")>1E-10 && nHitsPerChipValid[0]>20 && nHitsPerChipValid[1]>20",
				"goff",
				{x+"fitLast[0][1](100,-0.3,0.3)", x+"fitFirst[0][1](100,-0.3,0.3)", "telescopeFit(100,-0.3,0.3)"} )->Print((x+"partialfit.pdf(").c_str());
		combineHistogramsFromTree(*tree, { "meanDiffPerChipPerFitLast[1][0]."+x,	"meanDiffPerChipPerFitFirst[1][0]."+x, "meanDiffPerChip[1]."+x},
				"fabs(meanDiffPerChipPerFitLast[1][0]."+x+")>1E-10 && fabs(meanDiffPerChipPerFitFirst[1][0]."+x+")>1E-10 && nHitsPerChipValid[0]>20 && nHitsPerChipValid[1]>20",
				"goff",
				{x+"fitLast[1][0](100,-0.3,0.3)", x+"fitFirst[1][0](100,-0.3,0.3)", "telescopeFit(100,-0.3,0.3)"} )->Print((x+"partialfit.pdf").c_str());
		} catch(...) {};
		try{
		combineHistogramsFromTree(*tree, { "meanDiffPerChipPerFitLast[3][2]."+x,	"meanDiffPerChipPerFitFirst[3][2]."+x, "meanDiffPerChip[3]."+x},
				"fabs(meanDiffPerChipPerFitLast[3][2]."+x+")>1E-10 && fabs(meanDiffPerChipPerFitFirst[3][2]."+x+")>1E-10 && nHitsPerChipValid[2]>20 && nHitsPerChipValid[3]>20",
				"goff",
				{x+"fitLast[3][2](100,-0.3,0.3)", x+"fitFirst[3][2](100,-0.3,0.3)", "telescopeFit(100,-0.3,0.3)"} )->Print((x+"partialfit.pdf").c_str());
		combineHistogramsFromTree(*tree, { "meanDiffPerChipPerFitLast[2][3]."+x,	"meanDiffPerChipPerFitFirst[2][3]."+x, "meanDiffPerChip[2]."+x},
				"fabs(meanDiffPerChipPerFitLast[2][3]."+x+")>1E-10 && fabs(meanDiffPerChipPerFitFirst[2][3]."+x+")>1E-10 && nHitsPerChipValid[2]>20 && nHitsPerChipValid[3]>20",
				"goff",
				{x+"fitLast[2][3](100,-0.3,0.3)", x+"fitFirst[2][3](100,-0.3,0.3)", "telescopeFit(100,-0.3,0.3)"} )->Print((x+"partialfit.pdf)").c_str());
		} catch(...) {};
	}

}

void compareTelescopeAngles(std::string filename="test.root") {
	auto tree = getObjectFromFile<TTree>("fitTree", filename);
	std::string XZ{"XZ"};
	auto expressions = std::vector<std::string>{"(fitFirstThree."+XZ+".slope-fit."+XZ+".slope)", "(fitSecondThree."+XZ+".slope-fit."+XZ+".slope)", "fitFirstThree."+XZ+".slope-fitSecondThree."+XZ+".slope"};
	auto histnames = std::vector<std::string>{"firstAll(100,-0.002,0.002)", "lastAll(100,-0.002,0.002)","firstLast(100,-0.002,0.002)"};
//	auto titles = std::vector<std::string>{"(Difference of first and all)*#sqrt{2};Angle [Rad.]", "(Difference of last and all)*#sqrt{2}","Difference of first and last"};
	auto titles = std::vector<std::string>{"Difference of first and all;Angle [Rad.]", "Difference of last and all","Difference of first and last"};

	HistogramCombiner comb("telescopeAngles");
	for(int i=0; i<expressions.size(); i++) {
		auto h=getHistFromTree(*tree, expressions[i],"1", histnames[i], "goff");
		h->SetTitle(titles[i].c_str());
		h->SetTitle((h->GetTitle()+std::string(" (RMS = ")+to_string_precision(h->GetRMS()*1E3,3)+" mRad. )" ).c_str());
		comb.add( h );
	}

	comb.setNcolumn(1);
	comb.createCombined();
}
void scanTelescopeForFieldWires(std::string filename="test.root") {
	auto tree = getObjectFromFile<TTree>("fitTree", filename);
	auto canv=new TCanvas();
	canv->Print("scatteringInTelescope.pdf[");
	for(int i=150; i<230; i+=2) {
		tree->Draw(("fitFirstThree.XZ.slope-fitSecondThree.XZ.slope:fmod(fit.YZ.intercept+"+to_string(i)+"*fit.YZ.slope,2)").c_str(), "", "");
		canv->Update();
		canv->Print("scatteringInTelescope.pdf");
	}
	canv->Print("scatteringInTelescope.pdf]");


}

TProfile2D* drawPixelzCorrectionModulus(std::string filename="combinedFit.root") {
	auto th2 = getObjectFromFile<TProfile2D>("quad/zResidualByPixel", filename);
	int xbins=2, ybins=16;
	auto modulus=new TProfile2D{"modulus", "zResidualByPixel", xbins,0,double(xbins),ybins,0,double(ybins)};
	makeAppliedByPosition(th2, [&](double c, double x, double y){
				//std::cout<<x<<", "<<y<<": "<<c<<"\n";
				int ix(x), iy(y);
				if(fabs(c)<0.1)
					modulus->Fill(ix%xbins,iy%ybins,c);
				return c;
	});
	modulus->Draw("colz0");
	return modulus;
}

void compareOutOfSyncBackground() {
	CombineHistogramsFromFiles("eventhist/selectedHitAverageToTrackx",
			{"combinedFit_fiducial1Mnew.root","combinedFit_outOfSync1Mnew.root"},
			{"Selected tracks (offset=0)", "Telescope delayed (offset=1000 frames)"},
			[](HistogramCombiner &c){
				c.setNcolumn(1);
				c.setLogY();
			});

}

void drawOutOfSyncBackground(
	std::vector<std::string> files = {"combinedFit_inSync1Mnew.root","combinedFit_outOfSync1Mnew.root"} ) {

	std::vector<std::string> histnames{"tracks(60,-0.3,0.3)", "background(60,-0.3,0.3)"};
//	std::vector<std::string> histnames{"tracksPull(50,-5,5)", "background(50,-5,5)"};
	std::vector<std::string> titles {"Selected tracks", "Telescope delayed (offset=1000 frames)"};

	HistogramCombiner c{"outOfSync"};

	for(int i=0; i<files.size(); i++) {
//		auto h=getHistFromTreeInFile(files[i], "fitResults",
//				"meanQuadDiff.x/sqrt(pow(meanQuadError.x,2)+0.0288*0.0288)",
//				"Sum$(nHitsPerChipValid)<200 && Sum$(nHitsPerChipValid)>20 && matched&& fabs(centerDiffPerChipFit.x)>1E-10", histnames[i] ,"goff");
		auto h=getHistFromTreeInFile(files[i], "fitResults", "meanQuadDiff.x", "Sum$(nHitsPerChipValid)<200 && Sum$(nHitsPerChipValid)>20 && matched && fabs(centerDiffPerChipFit.x)>1E-10", histnames[i] ,"goff");
		h->SetTitle((titles[i]+";x-residual [mm];Events").c_str());
//		h->SetTitle((titles[i]+";Pull of x-residual;Events").c_str());
		c.add(h);
	}

	auto gausOffset=new TF1("gausOffset", "gaus(0)+[o]");
	gausOffset->SetParameters(1,0,1,1);
	c.hist[0]->Fit(gausOffset, "SM", "goff");
	gStyle->SetOptFit(0);

	c.setNcolumn(1);
	c.setLogY();
	c.setStyle(5);
	c.createCombinedOnCanvas(new TCanvas("outOfSync", "withOutOfSync", 600,600));
//	c.createCombined();

	auto stats=new StatsWrapper(0.15,0.73,0.4,0.68);
	stats->textSize=0.04;
//stats->add("#sigma", gausOffset->GetParameter(2), 2, "");
	stats->add("#sigma_{x}^{quad}", 1000*fabs(gausOffset->GetParameter(2)), 0, "#mum");
	stats->draw();
}

void drawQuadMeanResiduals( std::string fileName="combinedFit.root" ) {

	auto cuts="Sum$(nHitsPerChipValid)<200 && Sum$(nHitsPerChipValid)>20 && matched && fabs(centerDiffPerChipFit.x)>1E-10"; // && meanQuadPosition.z<8 && meanQuadPosition.z>4";
	auto gausOffset=new TF1("gausOffset", "gaus(0)+[o]");
	gausOffset->SetParameters(1,0,1,1);

	auto residual=getHistFromTreeInFile(fileName, "fitResults", "meanQuadDiff.x", cuts, "residual(60,-0.3,0.3)" ,"goff");
	residual->SetTitle(std::string(";x-residual [mm];Events").c_str());
	new TCanvas("meanQuadDiff.x");
	residual->Draw();
	residual->Fit(gausOffset);
//	HistogramCombiner residualHistogrammer{"residual", residual};
//	residualHistogrammer.makeLegend=false;
//	residualHistogrammer.setStyle(3);
//	residualHistogrammer.createCombined();


	auto pull=getHistFromTreeInFile(fileName, "fitResults", "meanQuadDiff.x/sqrt(pow(meanQuadError.x,2)+0.032*0.032)", cuts, "pull(50,-5,5)" ,"goff");
	pull->SetTitle(std::string(";x-residual/#sigma_{x};Events").c_str());
	new TCanvas("meanQuadDiff.x/meanQuadError.x");
	pull->Draw();
	pull->Fit(gausOffset);
//	HistogramCombiner pullHistogrammer{"pull", pull};
//	pullHistogrammer.makeLegend=false;
//	pullHistogrammer.setStyle(3);
//	pullHistogrammer.createCombined();

//
//	auto residualCenter=getHistFromTreeInFile(fileName, "fitResults", "centerDiffPerChipFit.x", cuts , "residualCenter(60,-0.3,0.3)" ,"goff");
//	residualCenter->SetTitle(std::string(";x-residual [mm];Events").c_str());
//	new TCanvas();
//	residualCenter->Draw();
//	gausOffset->SetParameters(1,0,1,1);
//	residualCenter->Fit(gausOffset);
//
//	auto pullCenter=getHistFromTreeInFile(fileName, "fitResults", "centerDiffPerChipFit.x/centerErrorPerChipFit.x", cuts, "pullCenter(50,-5,5)" ,"goff");
//	pullCenter->SetTitle(std::string(";x-residual/#sigma_{x};Events").c_str());
//	new TCanvas();
//	pullCenter->Draw();
//	pullCenter->Fit(gausOffset);

}

