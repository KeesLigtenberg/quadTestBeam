/*
 * BufferedTreeFiller.h
 *
 *  Created on: Nov 19, 2018
 *      Author: cligtenb
 */

#ifndef EVENTBUILDER_BUFFEREDTREEFILLER_H_
#define EVENTBUILDER_BUFFEREDTREEFILLER_H_

#include <vector>
#include <map>
#include <functional>
#include <iostream>

#include "TTree.h"

#include "Hit.h"

template <class K, class V>
class EntryBuffer {
public:
	void placeInBuffer(K tpcEntryNumber, const V&);
	bool isInBuffer(K tpcEntryNumber) const;
	V& getFromBuffer(K tpcEntryNumber);
	void removeFromBuffer(K tpcEntryNumber);
	void writeBufferUpTo(K tpcEntryNumber, std::function<void(V&)> function);
	void writeBuffer(std::function<void(V&)> function);
private:
	std::map<K,V> buffer; //buffer by tpcEntryNumber
};


template<class K, class V>
inline void EntryBuffer<K, V>::placeInBuffer(K tpcEntryNumber, const V& entry) {
	if(not buffer.insert( {tpcEntryNumber, entry} ).second) {
		std::cerr<<"Unexpected element in buffer!"; throw "Unexpected element in buffer!";
	}
}

template<class K, class V>
inline bool EntryBuffer<K, V>::isInBuffer(K tpcEntryNumber) const {
	if( buffer.find(tpcEntryNumber) == buffer.end() ) {
		return false;
	} else {
		return true;
	}
}

template<class K, class V>
inline V& EntryBuffer<K, V>::getFromBuffer(K tpcEntryNumber) {
	return buffer.at(tpcEntryNumber);
}

template <class K, class V>
inline void EntryBuffer<K,V>::removeFromBuffer(K tpcEntryNumber) {
	auto result=buffer.find(tpcEntryNumber);
	if(result==buffer.end()) { std::cerr<<"Could not find entry in buffer"; throw "Could not find entry in buffer"; }
	buffer.erase(result);
}

template <class K, class V>
inline void EntryBuffer<K,V>::writeBufferUpTo(K tpcEntryNumber, std::function<void(V&)> writeFunction) {
	for(auto e=buffer.begin();
			 e!=buffer.end() and e->first < tpcEntryNumber;
			 e=buffer.erase(e) ) {
		writeFunction(e->second);
	}
}

template <class K, class V>
void EntryBuffer<K,V>::writeBuffer(std::function<void(V&)> writeFunction) {
	for(auto& e : buffer) {
		writeFunction(e.second);
	}
	buffer.clear();
}

//helper class to fill tree buffered.
class BufferedTreeFiller {
public:
	BufferedTreeFiller() { setTreeBranches(); };

	struct TreeEntry {
		std::vector<Hit> chips[4]{};
		std::vector<int> nHitsShifted[4]{};
		long long triggerToA=0;
		unsigned triggerNumber=0;
	};

	void Write();
	int64_t GetEntries() { return tree.GetEntriesFast(); }
	void SetTreeDirectory(TFile* f) { tree.SetDirectory(f); }
	void placeInBuffer(int entryNumber, const TreeEntry& entry) {buffer.placeInBuffer(entryNumber, entry);}
	void emptyBufferUpTo(int entryNumber) { buffer.writeBufferUpTo(entryNumber, [this](TreeEntry&t){this->Fill(t);}); }
	void emptyBuffer() { buffer.writeBuffer([this](TreeEntry&t){this->Fill(t);}); }
	void removeFromBuffer(int entryNumber) { buffer.removeFromBuffer(entryNumber); }
	bool isInBuffer(int tpcEntryNumber) const { return buffer.isInBuffer(tpcEntryNumber); }
	TreeEntry& getEntry(int entryNumber) {
		if(not isInBuffer(entryNumber)) buffer.placeInBuffer(entryNumber, TreeEntry{});
		return buffer.getFromBuffer(entryNumber);
	};
	TTree& getTree() { return tree; }

private:
	TreeEntry currentEntry;
	EntryBuffer<int, TreeEntry> buffer;

	TTree tree{"data", "tree with run data"};

	void Fill(const TreeEntry& entry){
		currentEntry=entry;

		//count number of hits/chip/shift
		for(int i=0; i<4;i++) {
			for(int j=0; j<20; j++) {
				currentEntry.nHitsShifted[i].emplace_back( std::count_if(
					currentEntry.chips[i].begin(), currentEntry.chips[i].end(), [&j](const Hit& h){
						return h.nShiftedTrigger==j;
					}));
			}
		}

		tree.Fill();
	}
	void setTreeBranches();
};

inline void BufferedTreeFiller::setTreeBranches() {
	tree.Branch( "triggerToA", &currentEntry.triggerToA );
	tree.Branch( "triggerNumber", &currentEntry.triggerNumber );
	for(int i=0; i<4; i++) {
		tree.Branch( ("chip"+std::to_string(i)).c_str() , "std::vector<Hit>", currentEntry.chips+i);
		tree.Branch( ("nHitsShifted"+std::to_string(i)).c_str(), "std::vector<int>", currentEntry.nHitsShifted+i);
	}
}


#endif /* EVENTBUILDER_BUFFEREDTREEFILLER_H_ */
