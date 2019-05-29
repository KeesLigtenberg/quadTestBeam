//root macro


#include "/user/cligtenb/rootmacros/getObjectFromFile.h"
#include "/user/cligtenb/rootmacros/getHistFromTree.h"
#include "/user/cligtenb/rootmacros/AllCombiner.h"
#include "/user/cligtenb/rootmacros/makeGraphsFromFits.h"
#include "/user/cligtenb/rootmacros/histogramOperations.h"

struct gasFiles {
	std::vector<std::string> fileNames;
	std::vector<std::string> fileTitles;
	std::vector<double> voltages;
}

T2K {
	{	"run770_clusters.root",
		"run771_clusters.root",
		"run772_clusters.root",
		"run773_clusters.root",
		"run774_clusters.root",
		"run775_clusters.root",
		"run776_clusters.root",
		"run777_clusters.root" },
	{	"280 V", "290 V", "300 V", "310 V", "320 V", "330 V", "340 V", "350 V" },
	{	280,290,300,310,320,330,340,350 }
	//{ 53.5, 94,139,182,224,275,359,506 }
},

T2Kdropped {
	{	"run770_clusters_droph.root",
		"run771_clusters_droph.root",
		"run772_clusters_droph.root",
		"run773_clusters_droph.root",
		"run774_clusters_droph.root",
		"run775_clusters_droph.root",
		"run776_clusters_droph.root",
		"run777_clusters_droph.root" },
	{	"280 V", "290 V", "300 V", "310 V", "320 V", "330 V", "340 V", "350 V" },
	{	280,290,300,310,320,330,340,350 }
	//{ 53.5, 94,139,182,224,275,359,506 }
},

isobutane20 {
	{	"run779_clusters.root",
		"run780_clusters.root",
		"run781_clusters.root",
		"run782_clusters.root",
		"run783_clusters.root",
		"run784_clusters.root",
		"run785_clusters.root",
		"run786_clusters.root",
		"run787_clusters.root" },
	{	"340 V", "350 V", "360 V", "370 V", "380 V", "390 V", "400 V", "410 V", "420 V" },
	{340,350,360,370,380,390,400,410,420 }
//	{65.5,105,140,168,187,200,209,216,220}
},

isobutane20HT {
	{//	"run818_clusters.root",
//	"run819_clusters.root",
//	"run820_clusters.root",
//	"run821_clusters.root",
//	"run822_clusters.root",
		"run823_clusters.root",
		"run824_clusters.root",
		"run825_clusters.root",
		"run826_clusters.root",
		"run827_clusters.root" },
	{//"330 V", "340 V", "350 V", "360 V", "370 V", 
	"380 V", "390 V", "400 V", "410 V", "420 V" },
	{//	330, 340,350,360,370,
	380,390,400,410,420}
};

std::string cuts="meanCol>40 && meanCol<216 && meanRow>40 && meanRow<216 && stdDevRow>7 && stdDevCol>7";// && stdDevRow<13 && stdDevCol<13";

//move to macro file?
void annotateGraph( TGraph* g, std::vector<std::string> titles, double xOffset, double yOffset) {
	auto px=g->GetX();
	auto py=g->GetY();	
	for(int i=0; i<g->GetN(); i++) {
		auto text = new TText(px[i]+xOffset, py[i], titles[i].c_str() );
		text->SetTextSize(0.03);
		text->SetTextFont(42);
		text->Draw(); 
	}
}

TGraph* drawNHits(const gasFiles& gas, int chip=0) {
	std::vector<TH1*> hists;
	for(int i=0; i<gas.fileNames.size(); i++) {
		std::cout<<gas.fileNames[i]<<"\n";
		TTree* tree=getObjectFromFile<TTree>("clusters",gas.fileNames[i]);
		auto h=getHistFromTree(*tree, "nHits", cuts+"&& clusterChip=="+to_string(chip), gas.fileNames[i]+"(400,0,800)", "goff");
		h->SetTitle((gas.fileTitles[i]+";Number of hits;Clusters").c_str());
		if(h->GetEntries()) hists.push_back(h);
	}

	HistogramCombiner comb{"number of hits", hists};
	auto canv=comb.createCombined();
	canv->Print(("chip"+to_string(chip)+".pdf").c_str());

//	new TCanvas();
	auto graph=makeGraphFromHistExpression( hists, gas.voltages, [](TH1* hist, double) {
		hist->GetXaxis()->SetRangeUser( hist->GetMean()-1*hist->GetStdDev(), hist->GetMean()+1.5*hist->GetStdDev() );
		return std::make_pair( getMeanFromGausCoreFit(*hist) , hist->GetMean()/sqrt(hist->GetEntries()) ); //mean/sqrt(N) as error
	} );
//	graph->Draw("A*");
	//dumpGraphContent(graph);
	return graph;
}

void combineNumberOfHits(const gasFiles& gas) {	
	AllCombiner<TGraph> combiner{"nHitsGraphs"};
	for(int i=0; i<4; i++) {
		combiner.add( drawNHits(gas, i), "Chip "+to_string(i)+";Voltage [V];Most probable number of hits" );
	}
	auto canv=combiner.createCombined();
	canv->Print("nHits.pdf");
}

TGraph* drawToT(const gasFiles& gas, int chip=0) {
	std::vector<TH1*> hists;
	for(int i=0; i<gas.fileNames.size(); i++) {
		//new TCanvas(gas.fileTitles[i].c_str(),gas.fileTitles[i].c_str());
		std::cout<<gas.fileNames[i]<<"\n";
		TTree* tree=getObjectFromFile<TTree>("clusters",gas.fileNames[i]);
		auto h=getHistFromTree(*tree, "ToT", cuts+"&& chip=="+to_string(chip), gas.fileNames[i]+"(200,0,200)", "");
		h->SetTitle((gas.fileTitles[i]+";Time over threshold [25 ns];Hits").c_str());
		hists.push_back(h);
	}

	HistogramCombiner comb{"Time over threshold chip "+to_string(chip), hists};
	auto canvCombined=comb.createCombined();
	canvCombined->Print(("ToT_chip"+to_string(chip)+".pdf").c_str());

//	new TCanvas();
	auto graph=makeGraphFromHistExpression( hists, gas.voltages, [](TH1* hist, double) {
		//hist->GetXaxis()->SetRangeUser( hist->GetMean()-1*hist->GetStdDev(), hist->GetMean()+1.5*hist->GetStdDev() );
		return std::make_pair( hist->GetMean() , hist->GetMean()/sqrt(hist->GetEntries()) ); //mean/sqrt(N) as error
	} );
//	graph->Draw("A*");
	return graph;
}

void combineToTGraphs(const gasFiles& gas) {
	AllCombiner<TGraph> ToTGraphs("ToTGraphs");
	static TF1 expo("expo", "expo(0)");
	gStyle->SetOptFit(0);
	expo.SetLineWidth(1);
	for( int i=0; i<4; i++ ) {
		auto graph=drawToT(gas, i);
		graph->Fit(&expo, "N");
		ToTGraphs.add( graph,  "Chip "+to_string(i)+" (slope = "+to_string_precision(expo.GetParameter(1),3)+");Voltage [V];Mean ToT [25 ns]" );
	}
	auto canv=ToTGraphs.createCombined();
	canv->Print("ToT.pdf");
}

void combineVariousToTShapes(int iChip=0) {
	HistogramCombiner combined("combined");
	std::vector<std::string> files = { "run787_clusters.root", "run784_clusters.root", "run827_clusters.root"};
	std::vector<std::string> titles = { "420 V", "390 V", "420 V High Threshold"};
//	std::vector<std::string> files = { "run785_clusters.root", "run782_clusters.root", "run825_clusters.root"};
//	std::vector<std::string> titles = { "400 V", "370 V", "400 V High Threshold"};
	for(int i=0; i<files.size(); i++) {
		TTree* tree=getObjectFromFile<TTree>("clusters",files[i]);
		auto h=getHistFromTree(*tree, "ToT", cuts+"&& chip=="+to_string(iChip), files[i]+"(200,0,200)", "goff");
		combined.add(h, titles[i]+"; ToT [25 ns];Hits");
	}
	combined.normalise(20);
	combined.createCombined();
}

TGraph* graphToTNHits(const gasFiles& gas, int iChip) {
	auto tot=drawToT(gas,iChip);
	auto nHits=drawNHits(gas,iChip);
	if(tot->GetN()!=nHits->GetN()) cerr<<"number of graph points does not match!\n";
	return new TGraph(tot->GetN(), tot->GetY(), nHits->GetY());
}


void compareToTNHits(int iChip=0) {
	AllCombiner<TGraph> combined("combined");
	auto T2KGraph=graphToTNHits( T2K, iChip);
	combined.add( T2KGraph, "T2K;Mean ToT [25 ns];Most probable number of hits");
	combined.add( graphToTNHits( isobutane20, iChip), "Isobutane-18%");
	combined.add( graphToTNHits( isobutane20HT, iChip), "Isobutane-18% with high threshold");
	combined.setStyle(7);
	combined.createCombined();
	//annotateGraph(T2KGraph, T2K.fileTitles, 0, -10);
}

//returns number of control, isolated pairs, (control)
std::vector<double> compareToTPairs(std::string fileName) {
	HistogramCombiner combination("ToTPairComparison"+fileName);
	TTree* tree=getObjectFromFile<TTree>("clusters",fileName);
	std::vector<std::string> expr { "controlPairToTs", "isolatedPairToTs", "randomPairToTs"}; //, "randomIsolatedPairToTs" };
	std::vector<double> nPairs(expr.size());
	for( int i=0; i<expr.size(); i++ ) {
		for(auto iPair : {"first", "second"} ) {
			auto h=getHistFromTree(*tree, expr[i]+"."+iPair, cuts, expr[i]+iPair+"(100,0,200)", "");
			combination.add(h);
			h->SetTitle(expr[i].c_str());
			nPairs[i]=h->GetEntries();
		}
	}
	combination.setStyle(4);
	combination.skipNEntryLegend=2;
	combination.normalise();
	combination.createCombined();
	return nPairs;
}

void graphNumberOfPairs(const gasFiles& gas) {
	std::vector<double> nIsolated, nControl;
	for(int i=0; i<gas.fileNames.size(); i++) {
		const auto& file = gas.fileNames[i];
		auto nPairs = compareToTPairs(file);
		nControl.push_back(nPairs[0]);
		nIsolated.push_back(nPairs[1]);
		
		gPad->Print((to_string(int(gas.voltages[i]))+"V.pdf").c_str());
	}

	TGraph* graphIsolated = new TGraph( gas.fileNames.size(), gas.voltages.data(), nIsolated.data() );
	TGraph* graphControl = new TGraph( gas.fileNames.size(), gas.voltages.data(), nControl.data() );

	AllCombiner<TGraph> nPairsCombination("nPairsCombination");
	nPairsCombination.add( graphControl, "Number of control pairs;Voltage [V];Number of pairs");
	nPairsCombination.add( graphIsolated, "Number of isolated pairs");
	nPairsCombination.setYRangeRatio({1,2});
	nPairsCombination.createCombinedRatio();
}

TGraph* drawNNeigbours(const gasFiles& gas, int chip=0) {
	std::vector<TH1*> hists;
	for(int i=0; i<gas.fileNames.size(); i++) {
		std::cout<<gas.fileNames[i]<<"\n";
		TTree* tree=getObjectFromFile<TTree>("clusters",gas.fileNames[i]);
		auto h=getHistFromTree(*tree, "nNeighbours", cuts+"&& clusterChip=="+to_string(chip), gas.fileNames[i]+"(400,0,800)", "goff");
		h->SetTitle((gas.fileTitles[i]+";Number of neigbours;Clusters").c_str());
		if(h->GetEntries()) hists.push_back(h);
	}

	HistogramCombiner comb{"number of hits", hists};
	auto canv=comb.createCombined();
	canv->Print(("neigbours_chip"+to_string(chip)+".pdf").c_str());

	//new TCanvas();
	auto graph=makeGraphFromHistExpression( hists, gas.voltages, [](TH1* hist, double) {
		//hist->GetXaxis()->SetRangeUser( hist->GetMean()-1*hist->GetStdDev(), hist->GetMean()+1.5*hist->GetStdDev() );
		return std::make_pair( hist->GetMean() , hist->GetMean()/sqrt(hist->GetEntries()) ); //mean/sqrt(N) as error
	} );
	//graph->Draw("AP");
	return graph;
}

void combineNHitsNNeigbours(const gasFiles& gas) {	
	AllCombiner<TGraph> combiner{"nHitsGraphs"};
	for(int i=0; i<4; i++) {
		combiner.add( drawNHits(gas, i), "Hits chip "+to_string(i)+";Voltage [V];Number per event" );
		combiner.add( drawNNeigbours(gas, i), "Neighbours chip "+to_string(i) );
	}
	combiner.setStyle(4);
	combiner.compareTwoRatio=true;
	combiner.setYRangeRatio({0,1});
	auto canv=combiner.createCombinedRatio();
	canv->Print("nHitsNeigbhours.pdf");
}


