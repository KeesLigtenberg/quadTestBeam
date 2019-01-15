//root macro
#include <cstdlib>

#include "TPolyLine.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TH3.h"

#include "/user/cligtenb/rootmacros/getObjectFromFile.h"
#include "/user/cligtenb/rootmacros/getHistFromTree.h"
#include "/user/cligtenb/rootmacros/AllCombiner.h"
#include "/user/cligtenb/rootmacros/histogramOperations.h"
#include "/user/cligtenb/rootmacros/StatsWrapper.h"

//#include "../rootmacros/getObjectFromFile.h"
//#include "../rootmacros/getHistFromTree.h"
//#include "../rootmacros/AllCombiner.h"
//#include "../rootmacros/histogramOperations.h"
//#include "../rootmacros/StatsWrapper.h"


#include "Alignment/Alignment.h"

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

void combineHistogramsForChips(std::string histName="nHits", std::string fileName="combinedFit.root") {
	HistogramCombiner combination(histName+"combination");
	for(int i=0; i<4; i++) {
		auto hist=getObjectFromFile<TH1>(chipDirectories[i]+"/"+histName, fileName);
//		hist->Draw(); gPad->Update(); cin.get();
		combination.add( hist , "chip "+std::to_string(i+1) );
	}
	//combination.normalise();
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
void plotHitMap(std::string file="combinedFit.root", std::string alignFile="align.dat", std::string hitMapName="positionHitMap_local") {
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

void plotTimeWalkResiduals( std::string filename="fitted.root" ) {

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
	
	simple->SetParameters(0.366,0.066,0);
	means->Fit(simple, "", "", minToT,2.5);
	double offset=simple->GetParameter("offset");
//	means=shiftY(means, -simple->GetParameter("offset") );

	means->Draw();	
	simple->SetParameters(0.366,0.066,0);
	means->Fit(simple, "", "", minToT,2.5);


	means->GetYaxis()->SetTitle("Mean z-residual [mm]");
	means->GetYaxis()->SetTitleOffset(1.3);
	means->GetYaxis()->SetTitleSize(0.05);
	means->GetYaxis()->SetLabelSize(0.05);
	means->GetYaxis()->SetRangeUser(-1,4);
	means->GetXaxis()->SetTitleOffset(1.3);
	means->GetXaxis()->SetTitleSize(0.05);
	means->GetXaxis()->SetLabelSize(0.05);
	gPad->SetMargin(0.15,0.1,0.15,0.1);

	simple->SetLineWidth(1);
	gPad->SetTicks(1,1);
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	gStyle->SetOptTitle(0);

	auto stats=new StatsWrapper(0.5,0.6,0.88,0.88);
	stats->add("c_{1}", simple->GetParameter("c_{1}"), 3, "mm #mus");
	stats->add("t_{0}", simple->GetParameter("t_{0}"), 4, "#mus");
	stats->add("z_{offset}", offset, 3, "mm");
	//stats->addChiSquare(*simple);
	stats->draw();
	return means;
}
TH1* plotFittedTimeWalk( std::string filename="combinedFit.root", std::string object="quad/zResidualByToT") {
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
		comb.add(h, "chip "+to_string(i+1) );
	}
	comb.createCombined();
}


//fitSlicesY with specification of range
TH1* fitDiffusionSlices(TH2* h2, std::string x="z") {
//	TF1* gaus=new TF1("gaus","gaus(0)", -2,2);
	TF1* gausRange=new TF1("gausRange","gaus(0)", -0.6,0.6);
	gausRange->SetParameters(4E4,0.05,0.22);
//	TF1* exGaus=new TF1("exGaus", "[c]*[l]/2*exp([l]/2*(2*[m]+[l]*[s]*[s]-2*x))*TMath::Erfc( ([m]+[l]*[s]*[s]-x)/sqrt(2)/[s] )", -2, 2); //Exponentially modified gaussian distribution (see wiki)
//	exGaus->SetParameters(2E4,3.1,-0.25,0.2); // Constant, Lambda, Mean, Sigma
	//exGaus->FixParameter(1,3.1);

	h2->FitSlicesY(x=="z" ? gausRange : gausRange, 0/*firstbin*/, -1/*lastbin*/, 30/*min number of entries*/, "QNR");

	return dynamic_cast<TH1*>(gDirectory->Get( (h2->GetName()+std::string("_2")).c_str() ));
}


TF1* fitDiffusion( TH2* h2 , std::string x="x", double z0=-1, std::string canvname="canv") {
	double zmax=11;

	TF1* drift=new TF1("drift", ("sqrt( pow([#sigma_{"+x+"0}],2) + pow([D_{"+x+"}],2)/10*(x-[z0]) )").c_str(), z0, zmax);
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
	h2_2->GetYaxis()->SetRangeUser(0, x=="x"?0.5:0.6 );
	h2_2->GetXaxis()->SetRangeUser(z0,zmax);
	
	//guess parameters
	if(x=="z") drift->FixParameter(2,z0);
	else drift->SetParameter(2,z0); //z0
	drift->SetParameter(1,0.2); //D
	if(x=="x") 	drift->FixParameter(0, 0.0158771);
	else drift->SetParameter(0, 0.15);//sigma0

	//add error to fit
//	h2_2=addErrorToHist(h2_2, 1E-3); //set all error bins equal
//	h2_2->Fit(drift, "", "", z0, zmax);
	h2_2->Fit(drift, "", "", 1.5, zmax); //tmp until wiggle at low z is fixed!

	gStyle->SetOptTitle(0);	
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);

	h2_2->Draw();
	gPad->SetTicks(1,1);
	gPad->SetMargin(0.15,0.1,0.15,0.1);

	//nicer stat pane
	gStyle->SetOptFit(false);
	auto stats=new StatsWrapper();
	stats->add("D_{"+x+"}",  drift->GetParameter(1)*1E3, 0, "#mum/#sqrt{cm}" );
	if(x!="x") stats->add("#sigma_{"+x+"0}" ,  drift->GetParameter(0)*1E3, 0, "#mum" );
	stats->add("z0" ,  drift->GetParameter(2), 2, "mm" );
	stats->addChiSquare(*drift);
	stats->draw();
	
	//plot residuals
	const bool plotResiduals=false;
	if(plotResiduals) {
		auto residualHistogram=getResidualHistogram(h2_2, drift);
		new TCanvas(("residuals"+x).c_str(), "Residuals");
		residualHistogram->Draw();
	}
	
	return drift;

}

void plotDiffusionFromHist(std::string filename="combinedFit.root", std::string histogramName="ResidualByz_locExp") {
//	for(std::string chip : chipDirectories )
	std::string chip="quad";
	for(std::string x : {"x"}) {
		TH1* h=getObjectFromFile<TH1>( chip+"/"+x+histogramName, filename);
		TH2* h2=dynamic_cast<TH2*>(h);
		fitDiffusion(h2, x, -1, chip);
	}
}


void plotDiffusionCombined(std::string filename="fitted.root") {
	TH2D* histogram=nullptr;
	for(std::string x : {"x", "y", "z"}) {
		for(std::string chip : chipDirectories ) {
			TH2D* h=getObjectFromFile<TH2D>( chip+"/"+x+"ResidualByz", filename);
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


void plotSlicedDiffusionWithFit( std::string filename="combinedFit.root",  std::string object="quad/zResidualByzByToT_locExp" ) {
	auto hist3=getObjectFromFile<TH3D>(object, filename);
	auto zaxis=hist3->GetZaxis();

	double z0=-1.65;

	//HistogramCombiner slices("slices");
	std::vector<std::string> slicename = {"0.10 #mus < ToT < 0.30 #mus", "ToT > 0.30 #mus"};
	std::vector<std::pair<double,double> > binRanges = { {5,17}, {18, zaxis->GetNbins()+1} };
	std::vector<TH1*> histograms= { nullptr, nullptr };
	auto stats=new StatsWrapper();
	TF1* drift=new TF1("drift", "sqrt( pow([#sigma_{z0}],2) + pow([D],2)/10*(x-[z0]) )", 4.5, 23);;
	TLegend* legend = new TLegend();
	for(int i : {0,1} ) {
		zaxis->SetRange(binRanges[i].first, binRanges[i].second);
		auto proj=(TH2*)hist3->Project3D("yx");

		auto hsigma=fitDiffusionSlices(proj, "z");

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
			stats->add( "z0" ,  driftParams->GetParameter(2), 2, "mm" );

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
	legend->Draw();

	//slices.setStyle(10);slices.titleSize=0.05;
//	slices.setYRange({0,0.6});
//	slices.createCombined();

//	stats->draw();
	//gPad->SetMargin(0.15,0.1,0.15,0.1);

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



void fitDeformationCorrections(std::string fileName="combinedFit.root") {
	//	std::string dir="chip0";

	auto correction=new TF1("correction", "-[p1]/(1+pow((x-[p0])/[p2], 2)) + [p4]/(1+pow((x-[p3])/[p5],2))+ [p7]/(1+pow((x-[p6])/[p8],2))+[p9]");
	std::vector<std::vector<double> > estimatedValues={ {13.3,-2.05,1.067,13.48,-2.351,0.9421,-0.0166, 26, -2, 0.2} };

//	auto correction=new TF1("correction", "pol3");
//	std::vector<std::vector<double> > estimatedValues={ {1,1,1,1} };

//	struct { double min, max; } chipRange[4] = { {1.5,13}, {1.5,13}, {15.5,27}, {15.5,27} };
	struct { double min, max; } chipRange[4] = { {0,7}, {0,7}, {14,20}, {14,20} };
	for(int i=0; i<4; i++) {
		new TCanvas();
		auto xResidualByPosition=getObjectFromFile<TProfile2D>(chipDirectories[i]+"/xResidualByPosition_locExp", fileName );
		auto xResidualByx=xResidualByPosition->ProfileX();

		for(int j=0; j<correction->GetNpar(); j++) correction->SetParameter( j, estimatedValues[0][j] );
		xResidualByx->Fit(correction, "QS", "", chipRange[i].min, chipRange[i].max);

		std::cout<<"{";
		for(int j=0; j<correction->GetNpar()-1; j++) std::cout<<correction->GetParameter(j)<<", ";
		std::cout<<correction->GetParameter(correction->GetNpar()-1)<<"}\n";
	}
}

void compareDeformations(std::string fileName="combinedFit.root") {

	auto xResidualByPosition=getObjectFromFile<TProfile2D>("quad/xResidualByPosition_locExp", fileName );
	auto xResidualByx=xResidualByPosition->ProfileX();

	auto xResidualByPosition_cor=getObjectFromFile<TProfile2D>("quad/xResidualByPosition_corrected", fileName );
	auto xResidualByx_cor=xResidualByPosition_cor->ProfileX();

	HistogramCombiner deform("deform");
	deform.add(xResidualByx, "Before correction;x-position [mm]; x-residual [mm]");
	deform.add(xResidualByx_cor, "After correction");
	deform.createCombined();
}


//do frequency on an unweighted profile, i.e. all entreis should have same weight
TH1D* getFrequencyHistogramFromProfile(TProfile2D* original, double min=-0.1, double max=0.1, int nBins=80, double entryWeight=1.0) {
//	TH1D* frequencyHist=new TH1D(
//		(original->GetName()+std::string("freq")).c_str(),
//		(original->GetTitle()+std::string(" frequency;")+original->GetZaxis()->GetTitle()+";" ).c_str(),
//		nBins, min, max);
	TH1D* frequencyHistAll=new TH1D(
		(original->GetName()+std::string("freqAll")).c_str(),
		(original->GetTitle()+std::string(" frequency (with hits outside area);")+original->GetZaxis()->GetTitle()+";" ).c_str(),
		nBins, min, max);
//		TH1D* frequencyPull=new TH1D(
//			(original->GetName()+std::string("freqPull")).c_str(),
//			(original->GetTitle()+std::string(" frequencyPull;")+original->GetZaxis()->GetTitle()+";Frequency" ).c_str(),
//			100, -5, 5);
	for(int i=1;i<=original->GetNbinsX();i++) {
		for(int j=1; j<=original->GetNbinsY();j++) {
			int bin=original->GetBin(i,j);
			if( original->GetBinEntries(bin) <= 0 ) continue;
			double binContent= original->GetBinContent(i,j); //possible rounding error here
			double binError=original->GetBinError(i,j);
			frequencyHistAll->Fill(binContent);
//			if(!goodArea::isInsideArea( original->GetXaxis()->GetBinCenter(i), original->GetYaxis()->GetBinCenter(j) )) continue;
//			frequencyHist->Fill(binContent);
//				frequencyPull->Fill(binContent/error);

		}
	}
	return frequencyHistAll;
}

TProfile2D* getCombinedDeformations(std::string histName="quad/xResidualByPosition_locExp", std::vector<std::string> fileNames={"./run668/combinedFit.root","./run672/combinedFit.root","./run676/combinedFit.root"}) {
	auto prof=getObjectFromFile<TProfile2D>(histName, fileNames.front());
	for(int i=1; i<fileNames.size(); i++) {
		auto p=getObjectFromFile<TProfile2D>(histName, fileNames[i]);
		prof->Add(p);
	}

	return prof;
}

void combineDeformations(std::string histName="quad/xResidualByPosition_locExp", std::vector<std::string> fileNames={"./run668/combinedFit.root","./run672/combinedFit.root","./run676/combinedFit.root"}) {
	new TCanvas("hitmap", "hitmap", 900,1000)	;

	auto prof=getCombinedDeformations(histName, fileNames);

	//rebin
	prof->Rebin2D(4,4);
	removeBinsWithFewerEntries(prof, 1000);

	prof->SetTitle("");
	prof->GetYaxis()->SetTitleOffset(1.3);
	prof->GetXaxis()->SetTitleOffset(1.1);
	prof->GetZaxis()->SetTitleOffset(1.6);
	prof->Draw("colz0");

	gStyle->SetPalette(kRainBow);
	gPad->SetMargin(0.1,0.2,0.1,0.05);
	gPad->SetTicks(1,1);

	drawChipEdgesLocal("./run668/align.dat");

	new TCanvas();
	getFrequencyHistogramFromProfile(prof)->Draw();

}

void combineDeformationFrequencies(
		std::vector<std::string> histNames={"quad/xResidualByPosition_locExp",  "quad/xResidualByPosition_corrected" },
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

