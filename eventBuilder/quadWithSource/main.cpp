
#include "TSystem.h"

#include "processSource.h"

int main(int argc, const char* args[]) {
	if(argc!=3) {
		std::cout<<"usage: "<<args[0]<<" <infile> <outfile>\n";
		return 0;
	}

 	gROOT->ProcessLine(".L processSource.h+");

	processSource(args[1], args[2]);
}
