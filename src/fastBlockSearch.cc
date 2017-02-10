// #include <iostream>
// #include <vector>
// #include <fstream>
// #include <iomanip>
// #include <sstream>
// #include <string>
// #include <set>
#include <iostream>
#include <fstream>
#include <iomanip>

#define VERBOSE_SEEDS 1

#include "properties.hh"
#include "pp_profile.hh"
#include "pp_fastBlockSearcher.hh"

// #define SEEDSIZE 6
// #define SEED_COLUMNSCORE 0.8
#define MAX_SEQSIZE (unsigned(-1)/2)

//#define TESTLEN 1000000
#define MINCOUNT 5

using namespace std;

// set this with --cutoff=x
double cutoff = 0.7; //
   
// struct IndexComp {
//     IndexComp(const PP::Column& v) : vec(v) {}
//     const PP::Column& vec;
//     bool operator() (int a, int b) {
// 	return vec.Q(a) < vec.Q(b);
//     }
// };

// struct PrefilterType {
//     PrefilterType() : seedhitcount(0), lastOffset(-10) {}
//     int seedhitcount;
//     int lastOffset;
// };


void printBestResults(PP::Profile& prfl, PP::FsHitCollection& searchResults, double threshold, int offset) {
    multimap<double, PP::FsHitType> bestScores;
    searchResults.storeBestResults(bestScores, MINCOUNT, threshold);
    for (map<double,PP::FsHitType>::iterator it=bestScores.begin();
	 it != bestScores.end();
	 ++it) {
	vector<string> output(1, "--\n");
	cout << "Score:" << it->first << endl;
	for (PP::FsHitType* ht = &(it->second); ht != 0; ht = ht->predecessor) {
	    ostringstream sstrm;
	    sstrm << (ht->pos + offset) << "\t" << prfl[ht->blockNo].id;
	    if (ht->blockfrom != 0 || ht->blockto != prfl.blockSize(ht->blockNo)) 
		sstrm << "[" << ht->blockfrom << "," << ht->blockto << "]";
	    sstrm << "\t" << (ht->reverse ? '-' : '+') << "\t" << ht->score << "\t" 
		  << prfl[ht->blockNo].getBackDist().normed(ht->score) << "\n";
	    output.push_back(sstrm.str());
	}
	while (!output.empty()) {
	    cout << output.back();
	    output.pop_back();
	}
    }
}



int main(int argc, char* argv[]) {
    // filtering out blocks
    Properties::addProperty("/ProteinModel/blockpart_threshold_sens", "-1.0");
    Properties::addProperty("/ProteinModel/blockpart_threshold_spec", "8.0");
    Properties::addProperty("/ProteinModel/block_threshold_sens", "0.0");
    Properties::addProperty("/ProteinModel/block_threshold_spec", "6.0");
    bool abort = false;
    vector<string> argvec;
    for (int i=1; i<argc; i++) {
	string param(argv[i]);
	if (param.substr(0,2) == "--") {
	    int i = param.find("=");
	    if (param.find("=") != string::npos) {
		string key = param.substr(2, i-2);
		string value = param.substr(i+1);
		istringstream istrm(value);
		if (key == "cutoff") 
		    abort = !(istrm >> cutoff);
		else
		    Properties::addProperty(key, value);
	    } else
		abort = true;
	} else
	    argvec.push_back(param);
	if (abort) break;
    }

    if (abort || argvec.size() != 2) {
	cerr << "Usage: fastBlocksSearch [options] <seqs.fa> <prfl.ptn>\n";
	return -1;
    }

    // read in sequence file
//     map<string, SeqType> sequences;
    ifstream fstrm(argvec[0].c_str());
    if (!fstrm) {
	cerr << "Could not open \"" << argvec[0] << "\".\n";
	return -1;
    }
    
// Properties::addProperty("/ProteinModel/block_threshold_spec", "0.0");
    // Properties::addProperty("/ProteinModel/block_threshold_sens", "0.0");
    PP::initConstants();
    // PP::Profile::min_anchor_count = 0;

    // read in ptn
    PP::Profile prfl(argvec[1].c_str());
    
    int columnCount = 0;
    PP::FsSeedCollection seedColl(prfl);
    for (int b=0; b<prfl.blockCount(); b++)
	columnCount += prfl.blockSize(b);

    cerr << "Profile has " << prfl.blockCount() << " blocks with " << columnCount << " columns." << endl;
    cerr << "Total seed count: " << seedColl.size() << endl;
    cerr << "Reading sequence(s)";
    int predictionStart=0, predictionLen=MAX_SEQSIZE;
    try {
	predictionStart = Properties::getIntProperty("predictionStart") -1;
    } catch (...) {}
    try {
	predictionLen = Properties::getIntProperty("predictionEnd") - predictionStart;
    } catch (...) {}
    if (predictionLen <= 0) {
	cerr << "predictionStart must be less than predictionEnd!" << endl;
	exit(-1);
    }
    while (fstrm) {
	PP::FsHitCollection fastSearchResults(prfl.blockCount());

	string seqname, currSeq;
	getline(fstrm, seqname);
	if (seqname[0] != '>') {
	    cerr << "Not in FASTA format.\n";
	    return -1;
	}
	// int pos=seqname.find_first_of("\t\n\v\f\r ");
	// if (pos>=0)
	//     seqname.erase(pos);
	seqname = seqname.substr(1);
	
	int readlen = 0;
	int progress = 0;
	while (fstrm && fstrm.peek() != '>') {
	    string line;
	    getline(fstrm, line);
	    int pos = 0;
	    for (int i=0; i<line.length() && pos < predictionLen - currSeq.length(); i++) 
		if (isalpha(line[i])) //  || currSeq[i] == '-' || currSeq[i] == '.')
		    line[pos++] = tolower(line[i]);
	    line.erase(pos);
	    pos += readlen;
	    if (predictionStart <= pos) {
		if (predictionStart > readlen) 
		    currSeq += line.substr(predictionStart - readlen);
		else  
		    currSeq += line;
	    }
	    readlen = pos;
	    int newprogress = readlen/1000000;
	    while (newprogress > progress) {
		progress++;
		if (progress%10) 
		    cerr << ".";
		else
		    cerr << " " << progress << "M ";
	    }
	    if (currSeq.length() >= predictionLen) {
		if (currSeq.length() > predictionLen) 
		    cerr << "ERROR: prediction Length incorrect!" << endl;
		break;
	    }
	}
	cerr << endl << "Read " << readlen << " bps.";
	if (predictionLen < MAX_SEQSIZE) 
	    cerr << " Using " << predictionLen << ".";
	cerr << endl;

	PP::DNA::initSeq(currSeq);
	PP::CandidateCollection candidates(prfl, seedColl);
//	vector< deque<PrefilterType> >  prefilters[2];
// 	for (int b=0; b < prfl.blockCount(); b++) 
// 	    prefilters[0].push_back(deque<PrefilterType>(prfl[b].size() * 3, PrefilterType()));
// 	prefilters[1] = prefilters[0];
// 	vector<int> 
// 	    maxseedcount(prfl.blockCount(), 0), seedhitcounts(prfl.blockCount(), 0),
// 	    seedhitsizes(prfl.blockCount(), 0);
// 	vector<vector<int> > seedcountdist(prfl.blockCount(), vector<int>(6,0));
// 	deque<string> patterns(6, "XX");
	progress = 0;
	for (int t=2; t<currSeq.length(); t++) {
	    candidates.proceed(t, fastSearchResults);

	    int newprogress = 100 * t / currSeq.length();
	    while (newprogress > progress) {
		progress++;
		switch (progress%10) {
		    case 0 : cerr << progress; break;
		    case 9 : break;
		    default: cerr << ".";
		}
	    }
	}
	cerr << "100" << endl;
	cerr << " # hits: " << fastSearchResults.allHitCount() << endl; 
	cerr << " # hit groups: " << fastSearchResults.resultCount() << endl;

	cout << "Hits found in " << seqname << endl;
	printBestResults(prfl, fastSearchResults, cutoff * columnCount, predictionStart);
	cout << endl;
	if (predictionStart > 0 || predictionLen != MAX_SEQSIZE)
	    break;
    }
}
