#ifndef BINNEDCLUSTERER_H
#define BINNEDCLUSTERER_H

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

#include "HoughTransformer.h"
#include "../TrackFitter/PositionHit.h"
#include "../eventBuilder/Hit.h"
#include "DetectorConfiguration.h"

struct BinnedClusterer {

	constexpr static int nPlanes=6; //todo: make dynamic

	BinnedClusterer( double xmin, double xmax, double ymin, double ymax, int xbins, int ybins ) :
		xmin(xmin), xmax(xmax), ymin(ymin), ymax(ymax), xbins(xbins), ybins(ybins) {};

	virtual ~BinnedClusterer() {};

	typedef HoughTransformer::HitCluster HitCluster;

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

		//construct grid
		std::vector< std::vector< std::unique_ptr<HitCluster> > > houghGrid( xbins );
		for(auto& v : houghGrid) v.resize(ybins);

		constexpr bool DrawHistogram=false;
		std::unique_ptr<TH2D> graphicHistogram{
				DrawHistogram ?
				new TH2D("graphicHistogram", "Histogram of hough transform;x bin;y bin", xbins,0,xbins, ybins,0,ybins):
				nullptr };
		static TCanvas* canv=new TCanvas("houghTelCanv", "Canvas with cluster histogram", 600,400);

		for( unsigned plane=0; plane<hv.size(); plane++ ) {
			for(auto& h : hv[plane] ) {
				int binx= (h.position.x-xmin-angleOfTracksX*h.position.z)/(xmax-xmin)*xbins;
				int biny= (h.position.y-angleOfTracksY*h.position.z)/ymax*ybins;
				if( binx>=xbins ) binx=xbins-1;
				if( biny>=ybins ) biny=ybins-1;
				if( binx<0 ) binx=0;
				if( biny<0 ) biny=0;

				if(!houghGrid.at(binx).at(biny)) houghGrid.at(binx).at(biny)=std::unique_ptr<HitCluster>( new HitCluster() );
				houghGrid.at(binx).at(biny)->add(h, plane);

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

		//get grid positions and sort by size
		std::list< std::tuple<int, int, int> > gridPositions;
		for(int i=0; i<xbins; i++) {
			for(int j=0; j<ybins; j++) {
				if(houghGrid.at(i).at(j)) {
					int size=houghGrid.at(i).at(j)->clusterSize;
					if(size >= minCandidateSize) {
						gridPositions.emplace_back(size, i, j);
					}
				}
			}
		}
		gridPositions.sort();

		//merge neighbouring bins, starting with largest
		std::list<HitCluster> foundClusters;
		for(auto it=gridPositions.rbegin();	it!=gridPositions.rend(); it++) {
			int size, i, j;
			std::tie(size,i,j)=*it;

			if(!houghGrid.at(i).at(j)) continue;
			auto& currentCluster= *houghGrid.at(i).at(j);
			if(!currentCluster.clusterSize) continue;

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



#endif
