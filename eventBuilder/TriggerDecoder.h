//class to decode trigger numbers from spidr information
//Not programmed for efficiency

// use convention:
// example     __/-\___/----\___
// value:         1  2  4  8  16
// trigger bin:   0  1  2  3  4 
// rising edge:  0  1  2  3  4
//
// usage: construct TriggerDecoder with your number of bits and bin time length:
// binTimeWidth is bin of each timebin, firstBinOffset is the offset to the first rising edge
// TriggerDecoder td(15, 200, 900);
// Then call function with your rising edge times from the first rising edge in a trigger word
// auto myTriggerNumber= td.getNextTriggerFromTimes(myRisingEdgeTimes);
//
// Trigger number can only decrease in case trigger number 0, which is an empty list

#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <exception>

class TriggerDecoder {
public:
	TriggerDecoder(int nbits, double binTimeWidth=200, double firstBinOffset=950); // default values are for TLU used on June 8-9
	
	unsigned long long getNextTriggerFromTimes(const std::vector<int>& risingTimes);
    unsigned long long getNextTrigger(const std::vector<int>& risingedgepoints);

    //for exception
    struct TriggerDecoderException : std::exception {
    	TriggerDecoderException(std::string description) : exception(),  description(description) {}
    	const char* what() const noexcept {return description.c_str(); }
		std::string description;
    };

    //debug function
//    void setLastTrigger(unsigned short lt) { lastTrigger=lt; };


	double binTimeWidth;
	double firstBinOffset;


private:
	unsigned long long lastTrigger=0;
	const int maxTriggerIncrease=3;
	int triggersOutOfSync=0;
	const int maxTriggerOutOfSync=20;
	const int nbits;
	const double OutOfSyncIncreaseFactor=1.;
	
};

inline TriggerDecoder::TriggerDecoder(int nbits, double binTimeWidth, double firstBinOffset) :
		nbits(nbits) ,
		binTimeWidth(binTimeWidth),
		firstBinOffset(firstBinOffset)
{ 
	if(nbits>16) throw TriggerDecoderException("class currently has an implementation with maximum 16 bits!");
}


inline unsigned long long TriggerDecoder::getNextTriggerFromTimes(const std::vector<int>& risingTimes) {
	std::vector<int> risingEdgePoints;
	for(int time : risingTimes) {
		int point = std::floor( (time-firstBinOffset+binTimeWidth/2)/binTimeWidth );
		risingEdgePoints.push_back( point );
	}
	return getNextTrigger(risingEdgePoints);
}

//return true if combination is found, return false if reached end with wrong combination
inline bool findLowest(std::vector<char>& word, int pos, int lastPart ) {
	//check if done
	if(lastPart<=0) {
		//set all remaining bits to zero where allowed
		for(int i=int(word.size())-1; i>=0; i--) {
			if(!word[i]) {
				if(i==int(word.size())-1 or word[i+1]=='0') word[i]='0';
				else word[i]='1';
			}
		}
		for(char& c : word) if(!c) c='0';
		return true;
	}

	//reached end, but has not found the combination
	if(pos<0) return false;

	 //find first unknown position
	while(word[pos]) {
		//check if still in range
		if(--pos<0) return false;
	}

	//iterate 0,1 over unknow position starting with 0
	const int val=std::pow(2,pos);
	if(pos+1==int(word.size()) || word[pos+1]=='0' ) { //make sure not to introduce a new edge
		word[pos]='0';
		if(findLowest(word, pos-1,lastPart)) return true;
	}
	word[pos]='1';
	if(findLowest(word,pos-1, lastPart-val)) return true;

	//neither combination worked, reset and return
	word[pos]=0;
	return false;
}

std::string wordToString(std::vector<char> word) {
	std::string out;
	for(char c : word) {
		if(c) out+=c;
		else out+=' ';
	}
	return out;
}

inline unsigned long long TriggerDecoder::getNextTrigger(const std::vector<int>& risingedgepoints) {

	unsigned long long result=0;
	unsigned long long modulus=lastTrigger-lastTrigger%32768;
	if(!risingedgepoints.size()) {
		result=lastTrigger-lastTrigger%32768;
	} else {

		//put in information that is certain from rising edges into array
		std::vector<char> word(nbits,0);
		//always '0' before edge and '1' after
		for(auto& p : risingedgepoints) {
			if(p<0 or p>=nbits) {
				if(++triggersOutOfSync<maxTriggerOutOfSync) {
					return ++lastTrigger;
				} else {
					throw TriggerDecoderException("rising edge out of range! p="+std::to_string(p));
				}
			}
			if (p) word[p-1]='0';
			word[p]='1';
		}
		//all zeros before first edge
		if(!word[0]) for(int i=0; i<nbits; i++) {
			if(!word[i]) word[i]='0';
			else if(word[i]=='0') break;
			else throw TriggerDecoderException("unexpected character in word!");
		}
	
		//print certain bits
//		std::cout<<"certain bits: "<<wordToString(word)<<std::endl;

		//Now we find the correct combination:
		//use lastPart as the remaining value that has at least to be taken care of
		int lastPart=(lastTrigger+1)%32768;
		//decrease the remaining part by the certain information
		for(int i=0, n=1; i<nbits; n*=2, i++) {
			if(word[i]=='1') lastPart-=n;
		}

		//now try to find a possible next trigger word starting with the lowest combination
		if( findLowest(word,nbits-1,lastPart) ) {
			//print all bits
//			std::cout<<"deduced bits: "<<std::string(word.begin(), word.end() )<<std::endl; //found it!
		} else {
			throw TriggerDecoderException("Could not find next word: findLowest returned false! (Did the triggernumber decrease?) ");
			//in place of above exception the lastTrigger can also be reset and findLowest tried again!
		}

		//lastly convert word to number
		for(int i=0, n=1; i<nbits && word[i]; n*=2, i++) {
			if(word[i]=='1') result+=n;
			else if( !word[i] ) throw TriggerDecoderException("empty position in word!");
		}
//		std::cout<<"result="<<result<<" mod="<<modulus<<"\n";
		result=result%32768+modulus;
	}
	
	//do some checks
//	std::cout<<"result="<<result<<" last="<<lastTrigger<<"\n";
	if(result%32768 < lastTrigger%32768) {
		if(++triggersOutOfSync<maxTriggerOutOfSync) {
			if( int(lastTrigger%32768) <= 32768-maxTriggerIncrease) {
				if(!(result%32768) ) {
					throw TriggerDecoderException("wrong zero!");
				}
				std::cout<<"\nwarning: decrease in trigger number ignored: "<<result-lastTrigger<<"\n";
			}
			result=(lastTrigger+1);
		} else {
			throw TriggerDecoderException("trigger number decreased from "+std::to_string(lastTrigger)+" to "+std::to_string(result)+"!");
		}
	}
	else if( result==lastTrigger+1) triggersOutOfSync=0; //succesful!
	else if ( result<=lastTrigger+maxTriggerIncrease+OutOfSyncIncreaseFactor*triggersOutOfSync and triggersOutOfSync ) { std::cout<<"\nwarning: trigger increased by "<<result-lastTrigger<<"\n"; }
	else {
		if(lastTrigger==0) {//special case: this is the before the first actual trigger
			throw TriggerDecoderException("wrong zero!");
		} else if(++triggersOutOfSync<maxTriggerOutOfSync) {
			std::cout<<"\nwarning: increase in trigger number exceeded maxTriggerIncrease, ignoring increase "<<result-lastTrigger<<"\n";
			result=lastTrigger+1;
		} else {
			throw TriggerDecoderException("increase in trigger number exceeded maxTriggerIncrease and number of triggers out of sync, exceeds maxTriggerOutOfSync");
		}
	}
	
	//return result
	return lastTrigger=result;

}
