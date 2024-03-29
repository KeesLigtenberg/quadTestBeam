#ifndef GETHOUGHTRANSF_H
#define GETHOUGHTRANSF_H

#include <list>
#include <functional>
#include <memory>
#include <iostream>
#include <string>

#include "TStyle.h"
#include "TTree.h"
#include "TH2.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TView.h"
#include "TLegend.h"
#include "TPaletteAxis.h"
#include "TPaveLabel.h"
#include "TPaveStats.h"
#include "TSystem.h"
#include "TProfile2D.h"

#include "../TrackFitter/PositionHit.h"
#include "../eventBuilder/Hit.h"
#include "DetectorConfiguration.h"

struct HoughTransformer {

	constexpr static int nPlanes=6; //todo: make dynamic

	HoughTransformer( double xmin, double xmax, double ymin, double ymax, int xbins, int ybins ) :
		xmin(xmin), xmax(xmax), ymin(ymin), ymax(ymax), xbins(xbins), ybins(ybins) {};

	virtual ~HoughTransformer() {};

	struct HitCluster : std::list<PositionHit> { //todo: make member variable instead of inheritance (bad practice inheriting from stl classes)
		int clusterSize=0, planeHit[nPlanes] = {} ;
		void add( const PositionHit& h ) {
			push_back(h);
			++clusterSize;
			++planeHit[h.chip];
		};
		void add( const PositionHit& h, int plane ) {
			push_back(h);
			++clusterSize;
			++planeHit[plane];
		};
		void clear() {
			 std::list<PositionHit>::clear();
			 clusterSize=0;
			 for(int i=0; i<nPlanes; i++) planeHit[i]=0;
		};
		void mergeWith(HitCluster& o) {
			for(int i=0; i<nPlanes; i++) planeHit[i] += o.planeHit[i];
			clusterSize+=o.clusterSize;
			splice(begin(), o);
			o.clear();
		}
		int getNPlanesHit() const {
			int n=0;
			for(int i=0; i<nPlanes; i++) if(planeHit[i]) ++n;
			return n;
		}
		int getNHitsOnPlane(int i) const {
			return planeHit[i];
		}
		bool isPlaneHit(int i) const {
			return planeHit[i];
		}
		int recalculateNPlanesHit() {
			for(int i=0; i<nPlanes;i++) planeHit[i]=0;
			for( PositionHit& h : (*this) ) {
				++planeHit[h.chip];
			}
			return getNPlanesHit();
		}
		TVector3 getAveragePosition(bool rejectFlagged=false) const {
			TVector3 sum(0,0,0);
			int n=0;
			for(const auto& hit : *this) {
				if(rejectFlagged && hit.flag!=PositionHit::Flag::valid) continue;
				sum+=hit.position;
				n++;
			}
			sum*=(1./n);
			return sum;
		}
		int getNHitsUnflagged() {
			int n=0;
			for(const auto&h : (*this) ) {
				if(h.flag==PositionHit::Flag::valid) ++n;
			}
			return n;
		}
	};

	template<class T >
	static void drawClusters(const T& clusters, const DetectorConfiguration& detector);
	template<class T >
	static void drawCluster(const T& cluster, const DetectorConfiguration& detector);

	double xmin, xmax, ymin, ymax;
	int xbins, ybins;

	int minClusterSize=5, minCandidateSize=3;
	double angleOfTracksX=0., angleOfTracksY=0.;
	

	void changeStat(TH1&h) {
		gPad->Update();
		TPaveStats *st = (TPaveStats*)h.FindObject("stats");
		st->SetX1NDC(0.12); st->SetY1NDC(0.83);
		st->SetX2NDC(0.32); st->SetY2NDC(0.88);
		st->SetLineWidth(0);
		st->Draw();
		gPad->Update();
	}

	//telescope
	//vector of vectors, this means multiple planes
	virtual std::list<HitCluster> operator() ( const std::vector< std::vector<PositionHit> >& hv  ) {

		//list of pointers to all hits
		std::list<const PositionHit*> hitList;
		for(auto& v: hv) {
			for(auto& h: v) {
				hitList.push_back(&h);
			}
		}


		std::list<HitCluster> foundClusters;
		while(int(hitList.size())>minClusterSize) {

			//construct grid
			std::vector< std::vector< int > > houghGrid( xbins+1, std::vector<int>(ybins+1,0) );

			constexpr bool DrawHistogram=false;
			std::unique_ptr<TH2D> graphicHistogram{
					DrawHistogram ?
					new TH2D("graphicHistogram", "Histogram of hough transform;x bin;y bin", xbins,0,xbins, ybins,0,ybins):
					nullptr };
			static TCanvas* canv=DrawHistogram ? new TCanvas("houghTelCanv", "Canvas with cluster histogram", 600,400) : nullptr;

			for(auto& h : hitList ) {
				double angleRange=0.02;
				for(int binx=0; binx<=xbins; binx++) {
					double angle=angleOfTracksX-angleRange/2.+binx*angleRange/xbins;
					double r=h->position.x*cos(angle) + h->position.z*sin(angle);
					double rmin=0, rmax=20;
					int biny=(r-rmin)/(rmax-rmin)*ybins;

//					std::cout<<"fill r="<<r<<" and angle="<<angle<<"at bin "<<binx<<", "<<biny<<"\n";
					if(biny>=0 and biny<=ybins) houghGrid.at(binx).at(biny)++;

					if(DrawHistogram) graphicHistogram->Fill(binx,biny);
				}
			}

			//draw histogram
			if(DrawHistogram) {
				std::cout<<"drawing histogram of hough transform!"<<std::endl;
				canv->cd();
				gStyle->SetOptTitle(0);
				gStyle->SetOptStat("e");
				graphicHistogram->Draw("colz");
				changeStat(*graphicHistogram);
				canv->SetTicks(1,1);
				char c=std::cin.get();
				if(c=='w') {gPad->Print("telescopeCluster.pdf");}
				if(c=='q') { throw graphicHistogram;} //abuse of throw mechanism
			}

			//find largest cluster
			int maxBinEntries=0;
			int maxBinx=0, maxBiny=0;
			for(int col=0; col<=xbins; col++) {
				for(int row=0; row<=ybins; row++) {
					if(int(houghGrid[col][row]) > maxBinEntries) {
						maxBinEntries=houghGrid[col][row];
						maxBinx=col; maxBiny=row;
					}
				}
			}


			if(maxBinEntries<minClusterSize) break;

			//make houghcluster from largest
			HitCluster hc;
			for(auto ith=hitList.begin(); ith!=hitList.end(); ) {
				auto h=*ith;
				double angleRange=0.02;
				double angle=angleOfTracksX-angleRange/2.+maxBinx*angleRange/xbins;
				double r=h->position.x*cos(angle) + h->position.z*sin(angle);
				double rmin=0, rmax=20;
				int biny=(r-rmin)/(rmax-rmin)*ybins;
				if(biny==maxBiny) {
					hc.add(*h);
					ith=hitList.erase(ith);
				} else {
					ith++;
				}
			}
			foundClusters.push_back(hc);

		}

		return foundClusters;
	}
	

	//timepix
	//just a vector, this is one plane
	std::list<HoughTransformer::HitCluster> operator() ( const std::vector<PositionHit>& hv ) {

		//construct grid
		std::vector< std::vector< std::unique_ptr<HitCluster> > > houghGrid( xbins ); //Houghgrid[z][y]
		for(auto& v : houghGrid) v.resize(ybins);

		//set drawhistogram to true to make a histogram of the hough-like transform.
		constexpr bool DrawHistogram=false;
		static TCanvas* canv=nullptr;
		if(canv) canv->Clear();
		std::unique_ptr<TH2D> graphicHistogram{
				DrawHistogram ?
				new TH2D("graphicHistogram", "Histogram of hough transform;x bin;y bin", xbins,0,xbins, ybins,0,ybins ):
				nullptr};
		for(auto& h : hv ) {
			int binx= (h.position.x-xmin-angleOfTracksX*h.position.z)/(xmax-xmin)*xbins;
			int biny= (h.position.y-angleOfTracksY*h.position.z)/ymax*ybins;
			if( binx>=xbins ) binx=xbins-1;
			if( biny>=ybins ) biny=ybins-1;
			if( binx<0 ) binx=0;
			if( biny<0 ) biny=0;

			if(!houghGrid.at(binx).at(biny)) houghGrid.at(binx).at(biny)=std::unique_ptr<HitCluster>( new HitCluster() );
			houghGrid.at(binx).at(biny)->add(h);

			if(DrawHistogram) graphicHistogram->Fill(binx,biny);
		}

		//draw histogram
		if(DrawHistogram) {
			std::cout<<"drawing histogram of hough transform!"<<std::endl;
			if(!canv) canv=new TCanvas("houghCanv", "Canvas with cluster histogram", 600,400);
			canv->cd();
			graphicHistogram->Draw("colz");
			changeStat(*graphicHistogram);
			char c=std::cin.get();
			canv->SetTicks(1,1);
			gPad->Update();
			if(c=='w') {gPad->Print("timepixCluster.pdf");}
			if(c=='q') { throw 1;} //abuse of throw mechanism
		}


		struct ClusterCandidate {
			int size, x, y;
			ClusterCandidate(int size, int x, int y) : size(size),x(x),y(y) {};
		};

		//get grid positions
		std::list< ClusterCandidate > ClusterCandidatePositions;
		for(int i=0; i<xbins; i++) {
			for(int j=0; j<ybins; j++) {
				if(houghGrid.at(i).at(j)) {
					int size=houghGrid.at(i).at(j)->clusterSize;
					if(size >= minCandidateSize) {
						ClusterCandidatePositions.emplace_back(size, i, j);
					}
				}
			}
		}
		//sort by size
		ClusterCandidatePositions.sort([](const ClusterCandidate& a, const ClusterCandidate& b) {return a.size>b.size;});

		//merge neighbouring bins, starting with largest
		std::list<HitCluster> foundClusters;
		for(auto it=ClusterCandidatePositions.rbegin();	it!=ClusterCandidatePositions.rend(); it++) {
			int i=it->x, j=it->y;

			if(!houghGrid.at(i).at(j)) continue;
			auto& currentCluster= *houghGrid.at(i).at(j);
			if(!currentCluster.clusterSize) continue; //already merged with another

			//merge with neighbours
			for(int dx=-1; dx<=1; dx++) for(int dy=-1; dy<=1; dy++) {
				if(i+dx<0 or i+dx >= xbins or j+dy<0 or j+dy >= ybins or (!dx and !dy) ) continue; //outside grid
				if(! houghGrid.at(i+dx).at(j+dy) ) continue; //no hits in bin
				currentCluster.mergeWith( *houghGrid.at(i+dx).at(j+dy) );
			}

			if(currentCluster.clusterSize >= minClusterSize) {
				foundClusters.push_back( std::move(currentCluster) );
			}
			houghGrid.at(i).at(j).reset(nullptr);
		}

		return foundClusters;
	}




};

template<class T > //=std::list<HitCluster>
inline void HoughTransformer::drawClusters(const T& clusters, const DetectorConfiguration& detector) {
	static TCanvas* canvas=new TCanvas("clusters_canvas", "event display", 800,600);
	canvas->cd();
	TTree pointTree;
	PositionHit h;
	int cluster=1;
	pointTree.Branch("h", &h);
	pointTree.Branch("cluster", &cluster);

	for(auto& iClus : clusters) {
		for(auto& iHit : iClus ) {
			h=iHit;
			pointTree.Fill();
		}
		++cluster;
	}
	if(not pointTree.GetEntries()) return;

	pointTree.SetMarkerStyle(20);
	pointTree.Draw("h.position.y:h.position.z:h.position.x:cluster*10", "", "*colz");
//	pointTree.Draw("h.residual.y:h.position.z:h.residual.x:cluster*10", "", "*colz");


	gPad->Update();
	double theta=-20 /*70*/,phi=60;
//	std::cout<<"give angles!"<<std::endl;
//	std::cin>>theta>>phi;
	gPad->GetView()->RotateView(theta, phi);


	TH1* axisObject= dynamic_cast<TH1*>( gPad->GetPrimitive("htemp") );
	axisObject->GetXaxis()->SetLimits(detector.xmin(),detector.xmax());
	axisObject->GetZaxis()->SetLimits(detector.ymin(),detector.ymax());
//	axisObject->GetXaxis()->SetLimits(-.05,.05);
//	axisObject->GetZaxis()->SetLimits(-.05,.05);
	axisObject->DrawClone();
	gPad->Update();
}

template<class T > //=std::list<HitCluster>
inline void HoughTransformer::drawCluster(const T& cluster, const DetectorConfiguration& detector) {
	TTree pointTree;
	PositionHit h(1E4,1E4,1E4,0, Hit{0,0,0,0} );
	pointTree.Branch("h", "PositionHit", &h);
	pointTree.Fill(); //fill one with zero ToT to set scale

	for(auto& iHit : cluster ) {
		h=iHit;
//		std::cout<<int(h.ToT)<<"\n";
//		if(h.flag==PositionHit::Flag::shiftedTrigger) continue;
		if(h.flag==PositionHit::Flag::valid)
			pointTree.Fill();
	}
	if(not pointTree.GetEntries()) return;

	pointTree.SetMarkerStyle(20);
	gStyle->SetOptTitle(0);
	double totAxis=1.6;
	static TCanvas* canv=new TCanvas("cluster_display", "Event display", 800,600);
	canv->cd();
	pointTree.Draw( std::string("h.position.z:h.position.y:h.position.x:h.nShiftedTrigger").c_str() , "", "*colz"); //ToT to microseconds
//	pointTree.Draw( ("h.position.z:h.position.y:h.position.x:TMath::Min(h.ToT*0.025, "+std::to_string(totAxis)+")").c_str() , "", "*colz"); //ToT to microseconds
//	pointTree.Draw( std::string("h.position.z:h.residual.y:h.residual.x:h.chip").c_str() , "", "*colz"); //ToT to microseconds

	TH1* axisObject= dynamic_cast<TH1*>( gPad->GetPrimitive("htemp") );
	if(!axisObject) {std::cout<<"could not get axis object!?\n"; throw 1;}
	auto zaxis=axisObject->GetZaxis();
	if(!zaxis) {std::cout<<"could not get zaxis!?\n"; throw 1;}
	zaxis->SetLimits(detector.zmin(),detector.zmax());
	zaxis->SetTitle("z-axis (drift direction) [mm]") ;
	auto yaxis=axisObject->GetYaxis();
	if(!yaxis) {std::cout<<"could not get yaxis!?\n"; throw 1;}
	yaxis->SetTitle("y-axis (beam direction) [mm]");
	yaxis->SetLimits(detector.ymin(),detector.ymax());
	auto xaxis=axisObject->GetXaxis();
	if(!xaxis) {std::cout<<"could not get xaxis!?\n"; throw 1;}
	xaxis->SetLimits(detector.xmin(),detector.xmax());
	xaxis->SetTitle("x-axis [mm]");
	for(auto axis : {xaxis, yaxis} )
		axis->SetTitleOffset(1.3);
	for(auto axis : {xaxis, yaxis,zaxis} ) {
		axis->SetTitleSize(0.05);
		axis->SetLabelSize(0.05);
	}
//	axisObject->SetMaximum(totAxis);
//	axisObject->SetMinimum(0);

	axisObject->Draw("colz");
	gPad->SetMargin(0.1,0.175,0.15,0.1);//l r b t
	gPad->Update();

	TPaletteAxis* palette= dynamic_cast<TPaletteAxis*>(gPad->GetPrimitive("palette"));
	if(!palette) throw "could not find paletteAxis!";
	palette->SetX1NDC(0.85);
	palette->SetX2NDC(0.90);
	palette->SetY2NDC(0.74);
	palette->SetY1NDC(0.1);
	//draw TPaveText over Palette axis title
//	auto paletteAxisLabel = new TPaveLabel(0.96,0.1,1,0.75, "ToT [#mus]", "NDC");
	auto paletteAxisLabel = new TPaveLabel(0.96,0.1,1,0.75, "Frame [409.6 #mus]", "NDC");
	paletteAxisLabel->SetFillColor(kWhite);
	paletteAxisLabel->SetBorderSize(0);
	paletteAxisLabel->SetTextAngle(90);
	paletteAxisLabel->SetTextSize(0.08);
	paletteAxisLabel->SetTextFont(42);
	paletteAxisLabel->SetTextAlign(kHAlignCenter+kVAlignCenter);
	paletteAxisLabel->Draw();

	double theta=-20 /*-20*/,phi=60 /*10*/;
	//	std::cout<<"give angles!"<<std::endl;
	//	std::cin>>theta>>phi;
	gPad->GetView()->RotateView(theta, phi);

	TLegend* legend= new TLegend( 0.6, 0.8, 0.95,0.95 );
	legend->SetName("eventDisplayLegend");
	legend->AddEntry(axisObject, "Timepix hits", "p");
	axisObject->SetLineColor(kOrange+7);
	axisObject->SetLineWidth(2);
	legend->AddEntry(axisObject, "Telescope track", "l");
	legend->Draw();

	gPad->Update();

}

template<class T > //=std::list<HitCluster>
inline void drawCluster2D(const T& cluster, const DetectorConfiguration& detector) {
	TTree pointTree;
	PositionHit h(1E4,1E4,1E4,0, Hit{0,0,0,0} );
	pointTree.Branch("h", "PositionHit", &h);
//	pointTree.Fill(); //fill one with zero ToT to set scale

	for(auto& iHit : cluster ) {
		h=iHit;
//		std::cout<<int(h.ToT)<<"\n";
//		if(h.flag==PositionHit::Flag::shiftedTrigger) continue;
//		if(h.flag==PositionHit::Flag::valid)
			pointTree.Fill();
	}
	if(not pointTree.GetEntries()) return;

	pointTree.SetMarkerStyle(7);
	gStyle->SetOptTitle(0);
	double totAxis=1.6;
	static TCanvas* canv=new TCanvas("cluster_display2D", "Event display 2D", 600,800);
	canv->cd();
	pointTree.Draw( "h.position.y:h.position.x" , "", ""); //ToT to microseconds

	TH1* axisObject= dynamic_cast<TH1*>( gPad->GetPrimitive("htemp") );
	auto yaxis=axisObject->GetYaxis();
	if(!yaxis) {std::cout<<"could not get yaxis!?\n"; throw 1;}
	yaxis->SetTitle("y-axis (beam direction) [mm]");
	yaxis->SetLimits(detector.ymin(),detector.ymax());
	auto xaxis=axisObject->GetXaxis();
	if(!xaxis) {std::cout<<"could not get xaxis!?\n"; throw 1;}
	xaxis->SetLimits(detector.xmin(),detector.xmax());
	xaxis->SetTitle("x-axis [mm]");
	for(auto axis : {xaxis, yaxis} )
		axis->SetTitleOffset(1.1);
	for(auto axis : {xaxis, yaxis} ) {
		axis->SetTitleSize(0.05);
		axis->SetLabelSize(0.05);
	}
//	axisObject->SetMaximum(totAxis);
//	axisObject->SetMinimum(0);

	axisObject->Draw();
	pointTree.Draw( "h.position.y:h.position.x" , "", "same");
	gPad->SetMargin(0.15,0.1,0.1,0.05);//l r b t
	gPad->Update();

//
//	TLegend* legend= new TLegend( 0.6, 0.8, 0.95,0.95 );
//	legend->SetName("eventDisplayLegend");
//	legend->AddEntry(axisObject, "Timepix hits", "p");
//	axisObject->SetLineColor(kOrange+7);
//	axisObject->SetLineWidth(2);
//	legend->AddEntry(axisObject, "Telescope track", "l");
//	legend->Draw();

	gPad->Update();

}

template<class T > //=std::list<HitCluster>
inline void drawCluster2DPixel(const T& cluster) {
	TTree pointTree;
	PositionHit h(1E4,1E4,1E4,0, Hit{0,0,0,0} );
	pointTree.Branch("h", "PositionHit", &h);
//	pointTree.Fill(); //fill one with zero ToT to set scale

	TProfile2D prof("pixelProfile", "pixelProfile", 512,0,512,512,0,512);

	for(auto& iHit : cluster ) {
		h=iHit;
//		std::cout<<int(h.ToT)<<"\n";
		if(h.flag==PositionHit::Flag::shiftedTrigger) continue;
//		if(h.flag==PositionHit::Flag::valid)
//			pointTree.Fill();
				prof.Fill(
						(h.chip>=2)*256+(h.chip==0||h.chip==3)*(256-2*h.column)+h.column,
						(h.chip==1||h.chip==2)+255+(h.chip==0||h.chip==3)*-2*h.row+h.row,
						h.driftTime/4096.*25E-3);//choose between driftTime,
//						h.ToT*0.025); //and ToT
	}
//	if(not pointTree.GetEntries()) return;

	pointTree.SetMarkerStyle(7);
	gStyle->SetOptTitle(0);
	gStyle->SetPalette(kRainBow);
	double totAxis=5;
	static TCanvas* canv=new TCanvas("cluster_display2DPixel", "Event display 2D per pixel", 900,800);
	canv->cd();
//	pointTree.Draw( "h.ToT"
//			":(h.chip==1||h.chip==2)+255+(h.chip==0||h.chip==3)*-2*h.row+h.row"
//			":(h.chip>=2)*256+(h.chip==0||h.chip==3)*(256-2*h.column)+h.column"
//			">>hpixel(512,0,512,512,0,512)" , "", "profcolz");

	auto axisObject=&prof; //dynamic_cast<TH1*>( gPad->GetPrimitive("hpixel") );
	auto yaxis=axisObject->GetYaxis();
	if(!yaxis) {std::cout<<"could not get yaxis!?\n"; throw 1;}
	yaxis->SetTitle("Rows");
	auto xaxis=axisObject->GetXaxis();
	if(!xaxis) {std::cout<<"could not get xaxis!?\n"; throw 1;}
	xaxis->SetTitle("Columns");
	for(auto axis : {xaxis, yaxis} )
		axis->SetTitleOffset(1.1);
	for(auto axis : {xaxis, yaxis} ) {
		axis->SetTitleSize(0.045);
		axis->SetLabelSize(0.045);
	}
	axisObject->SetMaximum(totAxis);
	axisObject->SetMinimum(0);

	axisObject->DrawCopy("colz0");
	//pointTree.Draw( "h.position.y:h.position.x" , "", "same");
	gPad->SetMargin(0.15,0.15,0.15,0.05);//l r b t
	gPad->Update();

	gPad->Update();

}


//todo: move all draw functions away from houghtransformer
inline bool processDrawSignals() {
	static bool printedInfo=false;
	static bool pdfOpen=false;
	if(not printedInfo) { std::cout<<"<return> to continue, 'q' to break\n"; printedInfo=true; };

	char userSignal = std::cin.get();
//	std::cout<<"signal = "<<userSignal<<"\n";
	if (userSignal == 'q') {
		if(pdfOpen) {gPad->Print("eventDisplays.pdf]"); pdfOpen=false; }//close pdf
		return true;
	} else if (userSignal == 'l') {
		while (!gSystem->ProcessEvents()) {
			gSystem->Sleep(50);
		}
		return true;
	} else if (userSignal == 'w') {
		//rotate and write as animated gif!
		double phiView = 55;
		for (int thetaView = 0; thetaView < 360; thetaView += 2) {
			gPad->GetView()->RotateView(thetaView, phiView);
			gPad->Modified();
			gPad->Update();
			gPad->Print( thetaView == 358 ?	"eventAnimation.gif++5++" : "eventAnimation.gif+5");
//			gPad->Print(( "eventAnimation"+to_string(thetaView/2)+".png" ).c_str());
		}
	} else if (userSignal == 'a') {
		if(not pdfOpen) {
			gPad->Print("eventDisplays.pdf(");
			pdfOpen=true;
		} else {
			gPad->Print("eventDisplays.pdf");
		}
	}
	return false;
}

#endif
