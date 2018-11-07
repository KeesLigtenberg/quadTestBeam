#include "makeNoisyPixelMask.h"

#include <iostream>

#include "TH2.h"
#include "/user/cligtenb/rootmacros/getHistFromTree.h"
//#include "../../rootmacros/getHistFromTree.h"


#include "../TrackFitter/PositionHit.h"
//#include "getHistFromTree.h"

//mimosa
//threshold is maximum of number times mean
pixelMask makeNoisyPixelMask(TTree* hitTable, int plane, double threshold, std::pair<int,int> gridsize ) {
	int gridx=gridsize.first, gridy=gridsize.second;
	auto p=std::to_string(plane);
	TH2D hist("noisyPixelMaskHistogram", "histogram of all pixels",gridx,0,gridx,gridy,0,gridy );
	getHistFromTree(*hitTable, "mimosa["+p+"][].row:mimosa["+p+"][].column", "1", hist.GetName(),"goff", 1e5 /*DEBUG*/);

	auto mean = hist.GetEntries()/gridx/gridy;
//	std::cout<<"mean is "<<mean<<std::endl;
	pixelMask mask(gridx, std::vector<char>(gridy, 0) );
	int nmasked=0;
	for(int x=0; x<gridx; x++) {
		for(int y=0; y<gridy; y++) {
			if( hist.GetBinContent(x+1, y+1) > mean*threshold ) {
//				std::cout<<hist.GetBinContent(x+1, y+1)<<" > "<<mean*threshold<<std::endl;
				mask[x][y]=1;
				++nmasked;
			}
		}
	}

	std::cout<<"masked "<<nmasked<<" pixels for plane "<<plane<<std::endl;

	return mask;
}

//timepix
pixelMask makeNoisyPixelMask(TTree* hitTable, double threshold, std::pair<int,int> gridsize ) {
	int gridx=gridsize.first, gridy=gridsize.second;
	TH2D hist("noisyPixelMaskHistogram", "histogram of all pixels",gridx,0,gridx,gridy,0,gridy );
	getHistFromTree(*hitTable, "timepix[].row:timepix[].column", "1", hist.GetName() ,"goff",/* DEBUG*/ 1e5 );

	auto mean = hist.GetEntries()/gridx/gridy;
//	std::cout<<"mean is "<<mean<<std::endl;
	pixelMask mask(gridx, std::vector<char>(gridy, 0) );
	int nmasked=0;
	for(int x=0; x<gridx; x++) {
		for(int y=0; y<gridy; y++) {
			if( hist.GetBinContent(x+1, y+1) > mean*threshold ) {
//				std::cout<<hist.GetBinContent(x+1, y+1)<<" > "<<mean*threshold<<std::endl;
				mask[x][y]=1;
				++nmasked;
			}
		}
	}

	std::cout<<"masked "<<nmasked<<" pixels"<<std::endl;

	return mask;
}

template <class H=Hit>
std::vector<H> applyPixelMask(const pixelMask& mask, const std::vector<H>& hv ) {
	std::vector<H> newhv;
	for(const H& h:hv)	{
		if(mask.size()<=h.column) {std::cerr<<"applyPixelMask warning: column ("<<h.column<<") out of range ("<<mask.size()<<")\n"; continue;}
		if(mask[h.column].size()<=h.row) {std::cerr<<"applyPixelMask warning: row out of range\n"; continue;}
		if(! mask[h.column][h.row] ) {
			newhv.emplace_back(h);
		}
	}
//	std::cout<<"returning hv with size " << newhv.size()<<std::endl;
	return newhv;
}


//explicit instantations
template
std::vector<Hit> applyPixelMask(const pixelMask& mask, const std::vector<Hit>& hv );
template
std::vector<PositionHit> applyPixelMask(const pixelMask& mask, const std::vector<PositionHit>& hv );




