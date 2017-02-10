/*****************************************************************************\
 * Filename : pp_fastBlockSearcher.hh
 * Author   : Oliver Keller
 * Project  : Gene Prediction with Protein Family Patterns
 *
 * Description: A fast block search class to determine sequence parts
 *              relevant for the ppx extension
 *
 *
 * Date       |   Author              |  Changes
 *------------|-----------------------|---------------------------------
 * 09.09.09   | Oliver Keller         | creation of the file
 \******************************************************************************/

#ifndef __PP_BLOCKSEARCHER_HH
#define __PP_BLOCKSEARCHER_HH

// project includes
#include "pp_profile.hh"

// standard C/C++ includes
#include <deque>
#include <vector>
#include <map>
#include <ostream>
#include <iostream>
#include <iomanip>

using namespace std;

namespace PP {
    
    /*
     * FsSeedCollection
     */

    class FsSeedCollection {
	multimap<string, Position> seeds;
	
	static int expSeedCount;
	static double maxCoverage;
	multimap<string, Position>::iterator current;
	    
	static double blockP(const PP::Block& blk, string s, int from) {
	    double result = 1;
	    for (int i=0;  i<s.length(); i++) 
		result *= blk[from+i][GeneticCode::get_aa_from_symbol(s[i])];
	    return result;
	}
	static double blockQ(const Block& blk, string s, int from) {
	    double result = 1;
	    for (int i=0;  i<s.length(); i++) 
		result *= blk[from+i].Q(GeneticCode::get_aa_from_symbol(s[i])).doubleValue();
	    return result;
	}

    public:
	FsSeedCollection(const PP::Profile& prfl) {
	    for (int b = 0; b<prfl.blockCount(); b++) {
		const Block& blk = prfl[b];
#ifdef VERBOSE_SEEDS
		cerr << "Determining seeds for block " << b << " (" << blk.id << ", " << blk.size() << " columns).\n";
#endif
		int maxcount = expSeedCount / blk.size();

		for (int i=0; i<=blk.size()-3; i++) {
		    multimap<double,string> currtriples;
		    string s(3, 0);
		    for (int t1 = 0; t1<NUM_AA; t1++) {
			s[0] = GeneticCode::aa_symbols[t1];
			for (int t2 = 0; t2<NUM_AA; t2++) {
			    s[1] = GeneticCode::aa_symbols[t2];
			    for (int t3 = 0; t3<NUM_AA; t3++) {
				s[2] = GeneticCode::aa_symbols[t3];
				currtriples.insert(make_pair(blockQ(blk, s, i),s));
			    }
			}
		    }
		    int count = 0;
		    double p = 0.0;
		    map<double,string>::reverse_iterator it = currtriples.rbegin();
		    for (count = 0; count < maxcount && p < maxCoverage; count++) {
			p += blockP(blk, it->second, i);
			seeds.insert(make_pair(it->second, Position(b,i)));
			++it;
		    }
		}
	    }
	    current = seeds.end();
	}
	
	PP::Position* getFirst(string& s) {
	    current = seeds.lower_bound(s);
	    if (current == seeds.end())
		return 0;
	    if (current->first != s) {
		current = seeds.end();
		return 0;
	    }
	    return &current->second;
	}
	PP::Position* getNext() {
	    if (current == seeds.end())
		return 0;
	    const string& s = current->first;
	    if (++current == seeds.end())
		return 0;
	    if (current->first != s) {
		current = seeds.end();
		return 0;
	    }
	    return &current->second;
	}
	
	int size() {
	    return seeds.size();
	}
    }; // PP::FsSeedCollection

    /*
     * FsHitType
     */
    struct FsHitType;
    typedef deque<FsHitType*> HitQueue;

    struct FsHitType {
	FsHitType(int p, int b, bool r, PartScoreType sc) :
	    pos(p),
	    start(p),
	    blockNo(b),
	    reverse(r),
	    score(sc.score), 
	    blockfrom(sc.from),
	    blockto(sc.to),
	    pathScore(sc.score), 
	    predecessor(0) {}
	int pos;
	int start;
	int blockNo;
	bool reverse;
	double score;
	int blockfrom;
	int blockto;
	double pathScore;
	FsHitType* predecessor;
	void linkTo(HitQueue& queue) {
	    while(!queue.empty() && queue.front()->pos < pos - maxIntronLen) 
		queue.pop_front();
	    if (!queue.empty()) {
		predecessor = queue.front();
		start = predecessor->start;
		pathScore = predecessor->pathScore - intronMalus*(pos - predecessor->pos) + score;
	    }
	}
	void pushOn(HitQueue& queue) {
	    while (!queue.empty()) {
		FsHitType* ht = queue.back();
		if (ht->pathScore < pathScore + intronMalus*(pos - ht->pos))
		    queue.pop_back();
		else
		    break;
	    }
	    queue.push_back(this);
	}
	static int maxIntronLen;
	static double intronMalus;
    }; // PP::FsHitType

    /*
     * FsHitCollection
     */
    class FsHitCollection {
	vector<HitQueue> pendingHits[2];
	vector<FsHitType*> finalResult;
	vector<FsHitType*> allHits;
	int size;
    
    public:
	FsHitCollection(int n) : size(n) {
	    pendingHits[0].resize(n);
	    pendingHits[1].resize(n);
	}
	~FsHitCollection() {
	    for (int i=0; i<allHits.size(); i++) 
		delete allHits[i];
	}
	int allHitCount() {
	    return allHits.size();
	}
	int resultCount() {
	    return finalResult.size();
	}
	   
	void newHit(FsHitType*);
	void storeBestResults(multimap<double, FsHitType>&, int mincount, double threshold);
    }; // PP::FsHitCollection


    struct ScoringCandidate {
	ScoringCandidate() : inSeedCount(0), lastOffset(-10) {}
	int inSeedCount;
	int lastOffset;
    };

    struct CandidateCollection {
	typedef deque<ScoringCandidate> SeedCountsQueue;
	vector<SeedCountsQueue> theCandidates[2];
	FsSeedCollection& seedColl;
	const Profile& prfl;
	
	vector<int> maxseedcount; // maximum of inSeedCounts
	vector<int> seedhitsizes; // sum of inSeedCounts 
	vector<int> candcounts;   // for how many locations bestPartialLogScore is called
	vector<int> realhitcounts; // ...and successful
	
	// pending 2-mers of amino acids, one waiting for each frame/strand
	deque<string> patterns;


	CandidateCollection(const Profile& p, FsSeedCollection& coll) :
	    seedColl(coll),
	    prfl(p),
	    maxseedcount(p.blockCount(), 0),
	    seedhitsizes(p.blockCount(), 0),
	    candcounts(p.blockCount(), 0),
	    realhitcounts(p.blockCount(), 0),
	    patterns(6, "XX")
	{
	    for (int b=0; b<p.blockCount(); b++) 
		theCandidates[0].push_back(SeedCountsQueue(p.blockSize(b) * 3 -8, 
							ScoringCandidate()));
	    theCandidates[1] = theCandidates[0];
	}
	    
	void putSeedsToEntries(string currpat, bool reverse) {
	    for (PP::Position* it = seedColl.getFirst(currpat); it != 0; it = seedColl.getNext()) {
		PP::Position& pos = *it;
		int index = (reverse ? prfl.blockSize(pos.b) - pos.i - 3 : pos.i) * 3;
		    
		ScoringCandidate& pf = theCandidates[reverse][pos.b][index];
		if (pf.lastOffset >= 0 && (pf.lastOffset == pos.i || reverse == (pf.lastOffset < pos.i)))
		    cerr << "ERROR: last offset (" << pf.lastOffset << ") " << (reverse ? "<" : ">") 
			 << " current offset (" << pos.i << ") !" << endl;
		int range = abs(pf.lastOffset-pos.i);
		pf.inSeedCount += range < 3 ? range : 3;
		pf.lastOffset = pos.i;
	    }
	}

	void pushSeeds(int t) {
	    string nextPattern  
		= patterns.back() + GeneticCode::translate(PP::DNA::sequence + t-2);
	    patterns.pop_back();
	    string revPattern
		= GeneticCode::revtranslate(PP::DNA::sequence + t-2) + patterns.back();
	    patterns.pop_back();
	    patterns.push_front(nextPattern.substr(1));
	    patterns.push_front(revPattern.substr(0,2));
	    putSeedsToEntries(nextPattern, false);
	    putSeedsToEntries(revPattern, true);
	}

	void takeCandidates(int t, FsHitCollection& target, bool reverse) {
	    for (int b=0; b < prfl.blockCount(); b++) {
		theCandidates[reverse][b].push_front(ScoringCandidate());
		int hc = theCandidates[reverse][b].back().inSeedCount;
		theCandidates[reverse][b].pop_back();
		seedhitsizes[b] += hc;
		if (hc < 0) 
 		    cerr << "ERROR: seedhitcount < 0!\n";
		if (maxseedcount[b] < hc)
		    maxseedcount[b] = hc;
		if (hc > 4 + prfl.blockSize(b)/4) {
		    int blstart = t - prfl.blockSize(b) * 3;
		    if (blstart < 0)
			continue;
		    candcounts[b]++;
		    PP::PartScoreType pt;
		    prfl[b].bestPartialLogScore(reverse, blstart, pt);
		    if (pt.score >= 0) {
			FsHitType* newHit = new FsHitType(blstart, b, reverse, pt);
			target.newHit(newHit);
			realhitcounts[b]++;
		    }
		}
	    }
	}

	void proceed(int t, FsHitCollection& target) {
	    takeCandidates(t, target, false);
	    takeCandidates(t, target, true);
	    pushSeeds(t);
	}

	void debugging_output() {
	    for (int b=0; b<prfl.blockCount(); b++) {
		cerr << "Block " << b << " (" << prfl[b].id << ")" << endl;
		cerr << "    len: " << prfl.blockSize(b) << endl;
		cerr << "   av c: " << double(seedhitsizes[b])/DNA::len;
		cerr << "  max c: " << maxseedcount[b] << endl;
		cerr << "  % c>t: " << fixed << setw(4) << setprecision(2) << double(candcounts[b]*100)/DNA::len << endl;
	    }
	}
    }; // class PP::CandidateCollection
} // namespace PP
#endif // __PP_BLOCKSEARCHER_HH
