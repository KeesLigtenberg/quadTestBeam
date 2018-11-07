//root macro
#include <cstdlib>

#include "TPolyLine.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"

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
std::vector<std::string> chipDirectories={"chip1","chip2","chip3","chip4"};
std::vector<int> chipMapping={2,1,3,0};
std::vector<int> chipReverseMapping={4,2,1,3};//for canv->cd

void drawChipEdges(std::string alignFile="../align.dat") {
	Alignment	alignment(alignFile);
	for(int i=0; i<nChips; i++) {
		auto corners=alignment.getChipCorners(i);
		TPolyLine l;
		for(auto& corner : corners) {
			l.SetNextPoint(corner.x(), corner.y());
		}
		l.SetNextPoint(corners[0].x(), corners[0].y());
		l.DrawClone();
	}
	gStyle->SetOptStat(0);

}

void combineHistogramsForChips(std::string histName="nHits", std::string fileName="fitted.root") {
	HistogramCombiner combination(histName+"combination");
	for(auto& directory : chipDirectories) {
		auto hist=getObjectFromFile<TH1>(directory+"/"+histName, fileName);
//		hist->Draw(); gPad->Update(); cin.get();
		combination.add( hist , directory );
	}
	//combination.normalise();
	combination.createCombined();
}
void combineDrawHistogramsForChips(std::string expression, std::string cut, std::string drawOption, std::string histName="Hist", std::string fileName="fitted.root", std::string treeName="fitResults") {
	auto tree = getObjectFromFile<TTree>(treeName, fileName);
	HistogramCombiner combination(histName+"Combination");
	for(int i=0;i<4;i++) {
		auto hist=getHistFromTree(tree, expression, cut+" && chip=="+to_string(i), "chip"+to_string(i)+histName, drawOption.c_str());
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


void plotHitMap(std::string file="fitted.root", std::string alignFile="align.dat") {
	new TCanvas("hitmap", "hitmap", 850,1000)	;

	//hits
	auto first=true;
	for(const auto& dir:chipDirectories) {
		auto hist=getObjectFromFile<TH2>(dir+"/positionHitMap", file);
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
	drawChipEdges(alignFile);
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

	auto uncorrected=getObjectFromFile<TH2>("zResidualByToT", filename);
	TH1* uncorredtedHist=uncorrected->ProjectionY();
	auto corrected=getObjectFromFile<TH2>("zResidualByToTCorrected", filename);
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

void plotFittedTimeWalk(TH2* uncorrected) {
	uncorrected->FitSlicesY();
	auto means=dynamic_cast<TH1*>( gDirectory->Get("zResidualByToT_1") );
	if(!means) { cerr<<"failed to retrieve result from fitslicesy()\n"; return; };

	//rebin means
	{
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
	

	double minToT=0.15;
	auto simple=new TF1("2 parameters", "[c_{1}]/(x+[t_{0}])+[offset]", minToT,2.5);
	
	simple->SetParameters(0.366,0.066,0);
	means->Fit(simple, "", "", minToT,2.5);
	double offset=simple->GetParameter("offset");
	means=shiftY(means, -simple->GetParameter("offset") );

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
}
void plotFittedTimeWalk( std::string filename="fitted.root", std::string object="zResidualByToT") {
	auto uncorrected=getObjectFromFile<TH2>(object, filename);
	plotFittedTimeWalk(uncorrected);
}
//draw four timewalk on 4 pads in one canvas
void combineFittedTimeWalkForChips(std::string filename="fitted.root", std::string object="zResidualByToT") {
	auto canv=new TCanvas("canvas", "combined canvas", 1000,1000);
	canv->Divide(2,2);
	for(int i=0;i<4;i++) {
		canv->cd(i+1);
		plotFittedTimeWalk(filename, chipDirectories[i]+"/"+object);
	}
}


//fitSlicesY with specification of range
TH1* fitDiffusionSlices(TH2* h2, std::string x="z") {
	TF1* gaus=new TF1("gaus","gaus(0)", -2,2);
	TF1* gausRange=new TF1("gausRange","gaus(0)", -1,1);
	gausRange->SetParameters(4E4,0.05,0.22);
	TF1* exGaus=new TF1("exGaus", "[c]*[l]/2*exp([l]/2*(2*[m]+[l]*[s]*[s]-2*x))*TMath::Erfc( ([m]+[l]*[s]*[s]-x)/sqrt(2)/[s] )", -2, 2); //Exponentially modified gaussian distribution (see wiki)
	exGaus->SetParameters(2E4,3.1,-0.25,0.2); // Constant, Lambda, Mean, Sigma
	//exGaus->FixParameter(1,3.1);

	h2->FitSlicesY(x=="z" ? gausRange : gausRange, 0/*firstbin*/, -1/*lastbin*/, 5/*min number of entries*/, "QNR");

	return dynamic_cast<TH1*>(gDirectory->Get( (h2->GetName()+std::string("_2")).c_str() ));
}


TF1* fitDiffusion( TH2* h2 , std::string x="x", double z0=0, std::string canvname="canv") {
	TF1* drift=new TF1("drift", ("sqrt( pow([#sigma_{"+x+"0}],2) + pow([D_{"+x+"}],2)/10*(x-[z0]) )").c_str(), 4.5, 23);
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
	h2_2->GetYaxis()->SetRangeUser(0,1);
	h2_2->GetXaxis()->SetRangeUser(z0,25);
	
	//guess parameters
	drift->SetParameter(2,z0); //z0
	drift->SetParameter(1,0.2); //D
	drift->SetParameter(0, 0.15);//sigma0

	//add error to fit
	h2_2=addErrorToHist(h2_2, 1E-3); //set all error bins equal
	h2_2->Fit(drift, "", "", 0, 25);

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
	stats->add("#sigma_{"+x+"0}" ,  drift->GetParameter(0)*1E3, 0, "#mum" ); 
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

void plotDiffusionFromHist(std::string filename="fitted.root") {
	for(std::string chip : chipDirectories )
	for(std::string x : {"x", "y", "z"}) {
		TH1* h=getObjectFromFile<TH1>( chip+"/"+x+"ResidualByz", filename);
		TH2* h2=dynamic_cast<TH2*>(h);
		fitDiffusion(h2, x, 0, chip);
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

void plotZResidualsByToT(std::string filename="fitted.root") {
	HistogramCombiner combination("ToTCombination");
	TF1* gaus=new TF1("gaus","gaus(0)", -2,2);
	for(std::string chip : chipDirectories ) {
		TH2* h2=getObjectFromFile<TH2>( chip+"/zResidualByToT", filename);

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
	fitResults->Draw("XZ.intercept+XZ.slope*172:x");
	new TCanvas();
	fitResults->Draw("YZ.intercept+YZ.slope*172:z", "fabs(XZ.intercept+XZ.slope*172-x)<3");

}


