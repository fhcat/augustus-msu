/**********************************************************************
 * file:    exontrain.cc
 * licence: 
 *          
 * descr.:  training for exon model parameters
 * authors: Mario Stanke, mario@gobics.de
 *
 * date    |   author      |  changes
 * --------|---------------|------------------------------------------
 * 2.1.06  | Mario Stanke  | creating of file
 **********************************************************************/

#include "exonmodel.hh"

// project includes
#include "commontrain.hh"
#include "geneticcode.hh"
#include "gene.hh"
#include "projectio.hh"
#include "properties.hh"

// standard C/C++ includes
#include <fstream>
#include <sstream>


/*
 * ===[ ExonModel::buildModel ]===========================================
 */
void ExonModel::buildModel( const AnnoSequence* annoseq, int parIndex){
    if (!annoseq->anno) 
	throw ProjectError("Tried to train Exonmodel with 0 sequences.");
    
    if (verbosity)
	cout << " ** building model for exons *EXON*" << endl;
 
    gesbasen[0] = gesbasen[1] = gesbasen[2] = 0;
    for (int f=0; f<3; f++) {
	patterncount[f].assign(POWER4TOTHE( k+1 ), 0 );
    patternhighcount[f].assign(POWER4TOTHE( k+1 ), 0 );
    patternmediumcount[f].assign(POWER4TOTHE( k+1 ), 0 );
    patternlowcount[f].assign(POWER4TOTHE( k+1 ), 0 );
	initpatterncount[f].assign(POWER4TOTHE( k+1 ), 0 );
	highinitpatterncount[f].assign(POWER4TOTHE( k+1 ), 0 );
	mediuminitpatterncount[f].assign(POWER4TOTHE( k+1 ), 0 );
	lowinitpatterncount[f].assign(POWER4TOTHE( k+1 ), 0 );
	etpatterncount[f].assign(POWER4TOTHE( k+1 ), 0 );    
	highetpatterncount[f].assign(POWER4TOTHE( k+1 ), 0 );
	mediumetpatterncount[f].assign(POWER4TOTHE( k+1 ), 0 );
	lowetpatterncount[f].assign(POWER4TOTHE( k+1 ), 0 );
    }
    transInitMotif = new Motif(trans_init_window, tis_motif_memory, 1, tis_motif_radius);    // second parameter: k, war: *, 3, 1, 2
    if (etMotif){
	delete etMotif[0];
	delete etMotif[1];
	delete etMotif[2];
    }
    etMotif = new Motif*[3];
    etMotif[0] = new Motif(Constant::et_coding_len, etorder, etpseudocount);
    etMotif[1] = new Motif(Constant::et_coding_len, etorder, etpseudocount);
    etMotif[2] = new Motif(Constant::et_coding_len, etorder, etpseudocount);
//    delete baumwelch;
//    baumwelch = new BaumWelch(k, 3, numModels);

    ochrecount = ambercount = opalcount = 0;
    if (!hasLenDist) {
	lenCountSingle.assign(exonLenD+1,0);

	highlenCountSingle.assign(exonLenD+1,0);
	mediumlenCountSingle.assign(exonLenD+1,0);
	lowlenCountSingle.assign(exonLenD+1,0);

	lenCountInitial.assign(exonLenD+1,0);

	highlenCountInitial.assign(exonLenD+1,0);
	mediumlenCountInitial.assign(exonLenD+1,0);
	lowlenCountInitial.assign(exonLenD+1,0);

	lenCountInternal.assign(exonLenD+1,0);

	highlenCountInternal.assign(exonLenD+1,0);
	mediumlenCountInternal.assign(exonLenD+1,0);
	lowlenCountInternal.assign(exonLenD+1,0);

	lenCountTerminal.assign(exonLenD+1,0);

	highlenCountTerminal.assign(exonLenD+1,0);
	mediumlenCountTerminal.assign(exonLenD+1,0);
	lowlenCountTerminal.assign(exonLenD+1,0);

	numExonsOfType.assign(NUM_TYPES, 0);    
	numHugeExonsOfType.assign(NUM_TYPES, 0);    
    }
    const Gene* curgene;
    const AnnoSequence *as = annoseq;
    while (as){
      sequence = as->sequence;
      curgene = as->anno->genes; // assume here that sequences with multiple genes have already been split
      gweight = curgene->weight;
      if (curgene->clength % 3 == 0) {
	  try{
	      if (curgene->exons)
		  processExons( curgene );
	  }
	  catch( ExonModelError& exerr ){
	      cerr << "ExonModel::buildModel( " << curgene->id << " ):\t"
		   << exerr.getMessage() << endl << flush;
	  }
      } else {
	  if (verbosity)
	      cerr << "gene " << curgene->geneid << " transcr. " << curgene->id << " in sequence " << curgene->seqname << ": " 
		   << "coding length not a multiple of 3. Skipping..." << endl;
      }

      as = as->next;
    } 

/*
    // AA patterns
    if (!hasAAdep) {
	AADependency aasmall(1);
	aadep.setK(2);
	aadep.makeTransProbs(gene);
	aasmall.makeTransProbs(gene);
	 *
	 *  determine, whether a higher order markov chain brings different emission probabilities
	 *  by computing the Quotients P(A2|A1, A0) / P(A2|A1)
	 *
	for (int pn = 0; pn < aadep.trans.size(); pn++ ){
	    int shortpn;
	    shortpn = pn % ((int) (pow(20, aasmall.k+1) + 0.5));
	    aadep.trans[pn] = aadep.trans[pn] / aasmall.trans[shortpn];
	}
    }
*/

    buildProbabilities( );
    lastParIndex = parIndex;
}

/*
 * ===[ ExonModel::buildProbabilities ]===================================
 */ 
void ExonModel::buildProbabilities( ){
    //cerr << "entering buildProbabilities" << endl << flush;
 
    int numpatterns = patterncount[0].size();

 
    /*
     *  motif before translation initiation "ATG"
     */
    transInitMotif->makeProbs();
    //cout << "TransInit Motif" << endl;
    //transInitMotif->printProbs();
  
    /*
     * etMotif = motif before the donor splice site.
     * the index is the reading frame at the beginning (left) of the motif
     */
    etMotif[0]->makeProbs();
    etMotif[1]->makeProbs();
    etMotif[2]->makeProbs();    

   /*
    * Emissions probabilities and probabilities for shorter patterns are computed from frequencies 
    * of patterns of length k+1
    * Pls[i][f][p] is the probability of pattern number p of length (i+1) in frame f
    * The frame is the position of the last nucleotide. 
    */
    if (verbosity)
	cout << " number of bases in the reading frames: " 
	     << gesbasen[0]<< " " << gesbasen[1]<< " " << gesbasen[2]<< endl;
  
    Pls.assign( k+1, 3);
    highPls.assign( k+1, 3);
    mediumPls.assign( k+1, 3);
    lowPls.assign( k+1, 3);
 
    for (int f=0; f<3; f++) {
	Pls[k][f].resize(numpatterns);highPls[k][f].resize(numpatterns);mediumPls[k][f].resize(numpatterns);lowPls[k][f].resize(numpatterns);
	if (verbosity > 1)
	    cout << "--- frame = " << f << " ---    minPatSum = " << minPatSum << endl;
	if (minPatSum > 0) { 
	    determineShortPatterns(patterncount[f], k, minPatSum);
	    makeProbsFromCounts(Pls[k][f], patterncount[f], k, patpseudo, true);
	    determineShortPatterns(patternhighcount[f], k, minPatSum);
        makeProbsFromCounts(highPls[k][f], patternhighcount[f], k, patpseudo, true);
        determineShortPatterns(patternmediumcount[f], k, minPatSum);
        makeProbsFromCounts(mediumPls[k][f], patternmediumcount[f], k, patpseudo, true);
        determineShortPatterns(patternlowcount[f], k, minPatSum);
        makeProbsFromCounts(lowPls[k][f], patternlowcount[f], k, patpseudo, true);
	} else 
	    makeProbsFromCounts(Pls[k][f], patterncount[f], k, patpseudo, false);
        makeProbsFromCounts(highPls[k][f], patternhighcount[f], k, patpseudo, false);
        makeProbsFromCounts(mediumPls[k][f], patternmediumcount[f], k, patpseudo, false);
        makeProbsFromCounts(lowPls[k][f], patternlowcount[f], k, patpseudo, false);
    }
    
    // initial patterns, only emission probabilities
    for (int f=0; f<3; f++) {
	vector<Double> tempinitpatprobs(numpatterns);
	initemiprobs[f].resize(numpatterns);
	if (verbosity > 1)
	    cout << "--- initial frame = " << f << " ---    minPatSum = " << minPatSum << endl;
	if (minPatSum > 0) { 
	    determineShortPatterns(initpatterncount[f], k, minPatSum);
	    makeProbsFromCounts(tempinitpatprobs, initpatterncount[f], k, patpseudo, true);
	} else 
	    makeProbsFromCounts(tempinitpatprobs, initpatterncount[f], k, patpseudo, false);
	    
	    // make emission probabilities of the pattern probabilities and discard the latters
	computeEmiFromPat(tempinitpatprobs, initemiprobs[f], k);	
    }
    // end of initial patterns

    // initial patterns, only emission probabilities
    for (int f=0; f<3; f++) {
	vector<Double> hightempinitpatprobs(numpatterns);
	highinitemiprobs[f].resize(numpatterns);
	if (verbosity > 1)
	    cout << "--- initial frame = " << f << " ---    minPatSum = " << minPatSum << endl;
	if (minPatSum > 0) {
	    determineShortPatterns(highinitpatterncount[f], k, minPatSum);
	    makeProbsFromCounts(hightempinitpatprobs, highinitpatterncount[f], k, patpseudo, true);
	} else
	    makeProbsFromCounts(hightempinitpatprobs, highinitpatterncount[f], k, patpseudo, false);

	    // make emission probabilities of the pattern probabilities and discard the latters
	computeEmiFromPat(hightempinitpatprobs, highinitemiprobs[f], k);
    }
    // end of initial patterns

    // initial patterns, only emission probabilities
    for (int f=0; f<3; f++) {
	vector<Double> mediumtempinitpatprobs(numpatterns);
	mediuminitemiprobs[f].resize(numpatterns);
	if (verbosity > 1)
	    cout << "--- initial frame = " << f << " ---    minPatSum = " << minPatSum << endl;
	if (minPatSum > 0) {
	    determineShortPatterns(mediuminitpatterncount[f], k, minPatSum);
	    makeProbsFromCounts(mediumtempinitpatprobs, highinitpatterncount[f], k, patpseudo, true);
	} else
	    makeProbsFromCounts(mediumtempinitpatprobs, mediuminitpatterncount[f], k, patpseudo, false);

	    // make emission probabilities of the pattern probabilities and discard the latters
	computeEmiFromPat(mediumtempinitpatprobs, mediuminitemiprobs[f], k);
    }
    // end of initial patterns

    // initial patterns, only emission probabilities
    for (int f=0; f<3; f++) {
	vector<Double> lowtempinitpatprobs(numpatterns);
	lowinitemiprobs[f].resize(numpatterns);
	if (verbosity > 1)
	    cout << "--- initial frame = " << f << " ---    minPatSum = " << minPatSum << endl;
	if (minPatSum > 0) {
	    determineShortPatterns(lowinitpatterncount[f], k, minPatSum);
	    makeProbsFromCounts(lowtempinitpatprobs, lowinitpatterncount[f], k, patpseudo, true);
	} else
	    makeProbsFromCounts(lowtempinitpatprobs, lowinitpatterncount[f], k, patpseudo, false);

	    // make emission probabilities of the pattern probabilities and discard the latters
	computeEmiFromPat(lowtempinitpatprobs, lowinitemiprobs[f], k);
    }
    // end of initial patterns

    // internal exon terminal patterns, only emission probabilities
    for (int f=0; f<3; f++) {
	vector<Double> tempetpatprobs(numpatterns);
	etemiprobs[f].resize(numpatterns);
	if (verbosity > 1)
		cout << "--- internal exon terminal frame = " << f << " ---    minPatSum = " << minPatSum << endl;
	    if (minPatSum > 0) { 
		determineShortPatterns(etpatterncount[f], k, minPatSum);
		makeProbsFromCounts(tempetpatprobs, etpatterncount[f], k, patpseudo, true);
	    } else 
		makeProbsFromCounts(tempetpatprobs, etpatterncount[f], k, patpseudo, false);
	    // make emission probabilities of the pattern probabilities and discard the latters
	    computeEmiFromPat(tempetpatprobs, etemiprobs[f], k);	
    }
    // end of internal exon terminal patterns
    
    // internal exon terminal patterns, only emission probabilities
    for (int f=0; f<3; f++) {
	vector<Double> tempetpatprobs(numpatterns);
	highetemiprobs[f].resize(numpatterns);
	if (verbosity > 1)
		cout << "--- high internal exon terminal frame = " << f << " ---    minPatSum = " << minPatSum << endl;
	    if (minPatSum > 0) {
		determineShortPatterns(highetpatterncount[f], k, minPatSum);
		makeProbsFromCounts(tempetpatprobs, highetpatterncount[f], k, patpseudo, true);
	    } else
		makeProbsFromCounts(tempetpatprobs, highetpatterncount[f], k, patpseudo, false);
	    // make emission probabilities of the pattern probabilities and discard the latters
	    computeEmiFromPat(tempetpatprobs, highetemiprobs[f], k);
    }
    // end of internal exon terminal patterns

    // internal exon terminal patterns, only emission probabilities
    for (int f=0; f<3; f++) {
	vector<Double> tempetpatprobs(numpatterns);
	mediumetemiprobs[f].resize(numpatterns);
	if (verbosity > 1)
		cout << "--- internal exon terminal frame = " << f << " ---    minPatSum = " << minPatSum << endl;
	    if (minPatSum > 0) {
		determineShortPatterns(mediumetpatterncount[f], k, minPatSum);
		makeProbsFromCounts(tempetpatprobs, mediumetpatterncount[f], k, patpseudo, true);
	    } else
		makeProbsFromCounts(tempetpatprobs, mediumetpatterncount[f], k, patpseudo, false);
	    // make emission probabilities of the pattern probabilities and discard the latters
	    computeEmiFromPat(tempetpatprobs, mediumetemiprobs[f], k);
    }
    // end of internal exon terminal patterns

    // internal exon terminal patterns, only emission probabilities
    for (int f=0; f<3; f++) {
	vector<Double> tempetpatprobs(numpatterns);
	lowetemiprobs[f].resize(numpatterns);
	if (verbosity > 1)
		cout << "--- internal exon terminal frame = " << f << " ---    minPatSum = " << minPatSum << endl;
	    if (minPatSum > 0) {
		determineShortPatterns(lowetpatterncount[f], k, minPatSum);
		makeProbsFromCounts(tempetpatprobs, lowetpatterncount[f], k, patpseudo, true);
	    } else
		makeProbsFromCounts(tempetpatprobs, lowetpatterncount[f], k, patpseudo, false);
	    // make emission probabilities of the pattern probabilities and discard the latters
	    computeEmiFromPat(tempetpatprobs, lowetemiprobs[f], k);
    }
    // end of internal exon terminal patterns

    /* 
     * Compute lower order pattern probs from the ones above,
     * eg AAAA from AAAAA AAAAC AAAAG AAAAT
     */
    computeLowerOrderPats();

    emiprobs.probs[0].resize(Pls[k][0].size());
    emiprobs.probs[1].resize(Pls[k][1].size());
    emiprobs.probs[2].resize(Pls[k][2].size());
    highemiprobs.probs[0].resize(highPls[k][0].size());
    highemiprobs.probs[1].resize(highPls[k][1].size());
    highemiprobs.probs[2].resize(highPls[k][2].size());
    mediumemiprobs.probs[0].resize(mediumPls[k][0].size());
    mediumemiprobs.probs[1].resize(mediumPls[k][1].size());
    mediumemiprobs.probs[2].resize(mediumPls[k][2].size());    
    lowemiprobs.probs[0].resize(lowPls[k][0].size());
    lowemiprobs.probs[1].resize(lowPls[k][1].size());
    lowemiprobs.probs[2].resize(lowPls[k][2].size());
    
    emiprobs.order = k;
    computeEmiFromPat(Pls[k][0], emiprobs.probs[0], k);
    computeEmiFromPat(Pls[k][1], emiprobs.probs[1], k);
    computeEmiFromPat(Pls[k][2], emiprobs.probs[2], k);
    computeEmiFromPat(highPls[k][0], highemiprobs.probs[0], k);
    computeEmiFromPat(highPls[k][1], highemiprobs.probs[1], k);
    computeEmiFromPat(highPls[k][2], highemiprobs.probs[2], k);
    computeEmiFromPat(mediumPls[k][0], mediumemiprobs.probs[0], k);
    computeEmiFromPat(mediumPls[k][1], mediumemiprobs.probs[1], k);
    computeEmiFromPat(mediumPls[k][2], mediumemiprobs.probs[2], k);
    computeEmiFromPat(lowPls[k][0], lowemiprobs.probs[0], k);
    computeEmiFromPat(lowPls[k][1], lowemiprobs.probs[1], k);
    computeEmiFromPat(lowPls[k][2], lowemiprobs.probs[2], k);    

    if (!hasLenDist) 
	computeLengthDistributions();

    /*
     * Print out the numer of exons of each of the 8 types.
     * Only needed for making the transition matrix by hand.
     */
    if (verbosity)
	for (int type=singleGHigh; type <= terminalLow ; type++) {
	    cout << stateTypeNames[type] << " : " << numExonsOfType[type] << endl;
	}
    hasLenDist = true;
    if (verbosity>1){
	int totalcount = ochrecount + ambercount + opalcount;
	cout << "Frequency of stop codons:" << endl;
	if (totalcount>0) {
	    cout << "tag:" << setw(5) << ambercount << " (" << setprecision(3) << (float) ambercount/totalcount << ")" << endl;
	    cout << "taa:" << setw(5) << ochrecount << " (" << setprecision(3) << (float) ochrecount/totalcount << ")" << endl;
	    cout << "tga:" << setw(5) << opalcount << " (" << setprecision(3) << (float) opalcount/totalcount << ")" << endl;
	}
    }

    if (verbosity)
	cout << "end *EXON*" << endl;
}

/*
 * ===[ ExonModel::registerPars  ]===================================
 */
void ExonModel::registerPars (Parameters* parameters){
  Boolean CRFtrainCDS = true;
  try {
    CRFtrainCDS = Properties::getBoolProperty("CRFtrainCDS");
  } catch (...) {} 
  if (!parameters || !CRFtrainCDS)
    return;

  for (int idx=0; idx < Constant::decomp_num_steps; idx++)
    parameters->addMMGroup(&GCemiprobs[idx]);
}

/*
 * ===[ ExonModel::computeLowerOrderPats ]===================================
 */
void ExonModel::computeLowerOrderPats(){
   /* 
     * Compute lower order pattern probs Pls[i] (i=0..k-1) from the Pls[k]
     * eg AAAA from AAAAA CAAAA GAAAA TAAAA
     */
    int size;
    for( int i = k-1; i >= 0; i-- ){
	size = (int) (POWER4TOTHE( i+1 )+0.1);
        Pls[i][0].resize(size);
	Pls[i][1].resize(size);
	Pls[i][2].resize(size);

    highPls[i][0].resize(size);
    highPls[i][1].resize(size);
	highPls[i][2].resize(size);

    mediumPls[i][0].resize(size);
    mediumPls[i][1].resize(size);
	mediumPls[i][2].resize(size);

    lowPls[i][0].resize(size);
    lowPls[i][1].resize(size);
	lowPls[i][2].resize(size);

	for( int j = 0; j < size; j++ ){
            Double pi1 = Pls[i+1][0][j];
            Double pi2 = Pls[i+1][0][j + 1*size];
            Double pi3 = Pls[i+1][0][j + 2*size];
            Double pi4 = Pls[i+1][0][j + 3*size];
            Pls[i][0][j] = pi1+pi2+pi3+pi4;
            pi1 = Pls[i+1][1][j];
            pi2 = Pls[i+1][1][j + 1*size];
            pi3 = Pls[i+1][1][j + 2*size];
            pi4 = Pls[i+1][1][j + 3*size];
            Pls[i][1][j] = pi1+pi2+pi3+pi4;
	    pi1 = Pls[i+1][2][j];
            pi2 = Pls[i+1][2][j + 1*size];
            pi3 = Pls[i+1][2][j + 2*size];
            pi4 = Pls[i+1][2][j + 3*size];
            Pls[i][2][j] = pi1+pi2+pi3+pi4;
        }


	for( int j = 0; j < size; j++ ){
            Double pi1 = highPls[i+1][0][j];
            Double pi2 = highPls[i+1][0][j + 1*size];
            Double pi3 = highPls[i+1][0][j + 2*size];
            Double pi4 = highPls[i+1][0][j + 3*size];
            highPls[i][0][j] = pi1+pi2+pi3+pi4;
            pi1 = highPls[i+1][1][j];
            pi2 = highPls[i+1][1][j + 1*size];
            pi3 = highPls[i+1][1][j + 2*size];
            pi4 = highPls[i+1][1][j + 3*size];
            highPls[i][1][j] = pi1+pi2+pi3+pi4;
	    pi1 = highPls[i+1][2][j];
            pi2 = highPls[i+1][2][j + 1*size];
            pi3 = highPls[i+1][2][j + 2*size];
            pi4 = highPls[i+1][2][j + 3*size];
            highPls[i][2][j] = pi1+pi2+pi3+pi4;
    }

    for( int j = 0; j < size; j++ ){
			Double pi1 = mediumPls[i+1][0][j];
			Double pi2 = mediumPls[i+1][0][j + 1*size];
			Double pi3 = mediumPls[i+1][0][j + 2*size];
			Double pi4 = mediumPls[i+1][0][j + 3*size];
			mediumPls[i][0][j] = pi1+pi2+pi3+pi4;
			pi1 = mediumPls[i+1][1][j];
			pi2 = mediumPls[i+1][1][j + 1*size];
			pi3 = mediumPls[i+1][1][j + 2*size];
			pi4 = mediumPls[i+1][1][j + 3*size];
			mediumPls[i][1][j] = pi1+pi2+pi3+pi4;
		pi1 = mediumPls[i+1][2][j];
			pi2 = mediumPls[i+1][2][j + 1*size];
			pi3 = mediumPls[i+1][2][j + 2*size];
			pi4 = mediumPls[i+1][2][j + 3*size];
			mediumPls[i][2][j] = pi1+pi2+pi3+pi4;
            }

	for( int j = 0; j < size; j++ ){
		Double pi1 = lowPls[i+1][0][j];
		Double pi2 = lowPls[i+1][0][j + 1*size];
		Double pi3 = lowPls[i+1][0][j + 2*size];
		Double pi4 = lowPls[i+1][0][j + 3*size];
		lowPls[i][0][j] = pi1+pi2+pi3+pi4;
		pi1 = lowPls[i+1][1][j];
		pi2 = lowPls[i+1][1][j + 1*size];
		pi3 = lowPls[i+1][1][j + 2*size];
		pi4 = lowPls[i+1][1][j + 3*size];
		lowPls[i][1][j] = pi1+pi2+pi3+pi4;
	pi1 = lowPls[i+1][2][j];
		pi2 = lowPls[i+1][2][j + 1*size];
		pi3 = lowPls[i+1][2][j + 2*size];
		pi4 = lowPls[i+1][2][j + 3*size];
		lowPls[i][2][j] = pi1+pi2+pi3+pi4;
        }
    }
}


/*
 * ===[ ExonModel::storeGCPars ]=======================================
 * store GC dependent model probabilities in the respective array variables
 */
void ExonModel::storeGCPars(int idx){
  GCPls[idx] = Pls;
  
  highGCPls[idx] = highPls;
  mediumGCPls[idx] = mediumPls;
  lowGCPls[idx] = lowPls;
    
  GCemiprobs[idx] = emiprobs;
  
  highGCemiprobs[idx] = highemiprobs;
  mediumGCemiprobs[idx] = mediumemiprobs;
  lowGCemiprobs[idx] = lowemiprobs;
  
  std::stringstream out;
  out << "exon emiprob gc" << (idx+1);
  GCemiprobs[idx].setName(out.str());
  GCinitemiprobs[idx][0] = initemiprobs[0];
  GCinitemiprobs[idx][1] =  initemiprobs[1];
  GCinitemiprobs[idx][2] =  initemiprobs[2];

  highGCinitemiprobs[idx][0] = highinitemiprobs[0];
  highGCinitemiprobs[idx][1] =  highinitemiprobs[1];
  highGCinitemiprobs[idx][2] =  highinitemiprobs[2];

  mediumGCinitemiprobs[idx][0] = mediuminitemiprobs[0];
  mediumGCinitemiprobs[idx][1] =  mediuminitemiprobs[1];
  mediumGCinitemiprobs[idx][2] =  mediuminitemiprobs[2];

  lowGCinitemiprobs[idx][0] = lowinitemiprobs[0];
  lowGCinitemiprobs[idx][1] =  lowinitemiprobs[1];
  lowGCinitemiprobs[idx][2] =  lowinitemiprobs[2];


  GCetemiprobs[idx][0] =  etemiprobs[0];
  GCetemiprobs[idx][1] =  etemiprobs[1];
  GCetemiprobs[idx][2] =  etemiprobs[2];

  highGCetemiprobs[idx][0] =  highetemiprobs[0];
  highGCetemiprobs[idx][1] =  highetemiprobs[1];
  highGCetemiprobs[idx][2] =  highetemiprobs[2];

  mediumGCetemiprobs[idx][0] =  mediumetemiprobs[0];
  mediumGCetemiprobs[idx][1] =  mediumetemiprobs[1];
  mediumGCetemiprobs[idx][2] =  mediumetemiprobs[2];

  lowGCetemiprobs[idx][0] =  lowetemiprobs[0];
  lowGCetemiprobs[idx][1] =  lowetemiprobs[1];
  lowGCetemiprobs[idx][2] =  lowetemiprobs[2];

  if (transInitMotif)
    GCtransInitMotif[idx] =  *transInitMotif;
  if (etMotif){
    *GCetMotif[idx][0] = *etMotif[0];
    *GCetMotif[idx][1] = *etMotif[1];
    *GCetMotif[idx][2] = *etMotif[2];
  }
}

/*
 * printProbabilities
 * idx - the index of the GC content class
 * suffix - an identifyer to append to the parameter file name in case there are multiple variants generated during training
 */
void ExonModel::printProbabilities(int idx, BaseCount *bc, const char* suffix){
    Seq2Int s2i_e(k+1);

        try{ 
	    const char* fname = Properties::getProperty( "/ExonModel/outfile" );
	    string filename = Constant::fullSpeciesPath() + fname;
	    if (suffix)
		filename += suffix;
	    if (idx < 0) { // do nothing but deleting the outfile
		ofstream ofstrm(filename.c_str(), ios::trunc);
		ofstrm.close();
		return;
	    }
            ofstream ofstrm(filename.c_str(), ios::app);
	    if (verbosity)
  	        cout << "Writing exon model parameters [" << idx+1 << "] to file " << filename << "." << endl;
//	    LLDouble::setOutputPrecision(3);
	    if (idx == 0 ) { 
		ofstrm << "#exon model parameters\n# begin of content independent part" << endl;
		
		ofstrm << "\n# Length distributions\n[LENGTH]" << endl;
		ofstrm << "# maximal individually stored length probability =\n" << exonLenD << endl;
		ofstrm << "# slope of smoothing bandwidth =\n" << slope_of_bandwidth << endl;
		ofstrm << "# smoothing minwindowcount =\n" << minwindowcount << endl;
		ofstrm << "# length single  initial  internal  terminal" << endl;
		ofstrm << "# total number of exons of above types" << endl;

		ofstrm << "       " << numSingle << setw(15) << numInitial << setw(15)
				       << numInternal <<setw(15) << numTerminal << endl;


		ofstrm << "       " << highnumSingle << setw(15) << highnumInitial << setw(15)
		       << highnumInternal <<setw(15) << highnumTerminal << endl;

		ofstrm << "       " << mediumnumSingle << setw(15) << mediumnumInitial << setw(15)
		       << mediumnumInternal <<setw(15) << mediumnumTerminal << endl;

		ofstrm << "       " << lownumSingle << setw(15) << lownumInitial << setw(15)
		       << lownumInternal <<setw(15) << lownumTerminal << endl;


		ofstrm << "# number of exons exceeding length d" << endl;

		ofstrm << "       " << numHugeSingle << setw(15) << numHugeInitial << setw(15) 
		       << numHugeInternal << setw(15) << numHugeTerminal << endl;

		ofstrm << "       " << highnumHugeSingle << setw(15) << highnumHugeInitial << setw(15)
		       << highnumHugeInternal << setw(15) << highnumHugeTerminal << endl;

		ofstrm << "       " << mediumnumHugeSingle << setw(15) << mediumnumHugeInitial << setw(15)
		       << mediumnumHugeInternal << setw(15) << mediumnumHugeTerminal << endl;

		ofstrm << "       " << lownumHugeSingle << setw(15) << lownumHugeInitial << setw(15)
		       << lownumHugeInternal << setw(15) << lownumHugeTerminal << endl;

		
	/*	ofstrm << "# 1000 P(len=k), k=0,1,..., " << exonLenD << endl;
		for( int i = 0; i <= exonLenD; i++ ){
		    ofstrm << i 
			   << "\t" << 1000 * lenDistSingle[i] << "\t" << 1000 * lenDistInitial[i] 
			   << "\t" << 1000 * lenDistInternal[i] << "\t" << 1000 * lenDistTerminal[i] 
			   << endl;
		}*/

		ofstrm << "# 1000 P_high(len=k), k=0,1,..., " << exonLenD << endl;
		for( int i = 0; i <= exonLenD; i++ ){
		    ofstrm << i
			   << "\t" << 1000 * highlenDistSingle[i] << "\t" << 1000 * highlenDistInitial[i]
			   << "\t" << 1000 * highlenDistInternal[i] << "\t" << 1000 * highlenDistTerminal[i]
			   << endl;
		}

		ofstrm << "# 1000 P_medium(len=k), k=0,1,..., " << exonLenD << endl;
		for( int i = 0; i <= exonLenD; i++ ){
		    ofstrm << i
			   << "\t" << 1000 * mediumlenDistSingle[i] << "\t" << 1000 * mediumlenDistInitial[i]
			   << "\t" << 1000 * mediumlenDistInternal[i] << "\t" << 1000 * mediumlenDistTerminal[i]
			   << endl;
		}

		ofstrm << "# 1000 P_low(len=k), k=0,1,..., " << exonLenD << endl;
		for( int i = 0; i <= exonLenD; i++ ){
		    ofstrm << i
			   << "\t" << 1000 * lowlenDistSingle[i] << "\t" << 1000 * lowlenDistInitial[i]
			   << "\t" << 1000 * lowlenDistInternal[i] << "\t" << 1000 * lowlenDistTerminal[i]
			   << endl;
		}

/*
		ofstrm << "\n# actual lengths " << endl;
		for( int i = 0; i <= exonLenD; i++ ){
		    ofstrm << i 
			   << "\t" << lenCountSingle[i] << "\t" << lenCountInitial[i] 
			   << "\t" << lenCountInternal[i] << "\t" << lenCountTerminal[i] 
			   << endl;
		}

*/
		ofstrm << "# end of content independent part" << endl;
	    }

	    ofstrm << "\n# data set number\n[" << idx+1 << "]" << endl;
	    if (bc)
		ofstrm << "# (a,c,g,t)= " << *bc << endl;
            /*ofstrm << "#\n# Probabilities file for the exon model\n#\n"
                   << endl;


            ofstrm << "\n# Die P_l's\n[P_ls]" << endl;
            ofstrm << "# k = " << k << endl;
            for( int i = 0; i <= k; i++ ){
                ofstrm << "# l=\n" <<  i << endl;
                Seq2Int s2i(i+1);
                ofstrm << "# Values" << endl;
                for( int j = 0; j < GCPls[idx][i][0].size(); j++ )
                    ofstrm << s2i.inv(j) << "\t" << GCPls[idx][i][0][j]
                           << "\t     " << GCPls[idx][i][1][j] << "\t     "
                           << GCPls[idx][i][2][j] << endl;
            }*/

                ofstrm << "#\n# High Probabilities file for the exon model\n#\n"
                   << endl;


            ofstrm << "\n# High Die P_l's\n[High P_ls]" << endl;
            ofstrm << "# k = " << k << endl;
            for( int i = 0; i <= k; i++ ){
                ofstrm << "# l=\n" <<  i << endl;
                Seq2Int s2i(i+1);
                ofstrm << "# Values" << endl;
                for( int j = 0; j < highGCPls[idx][i][0].size(); j++ )
                    ofstrm << s2i.inv(j) << "\t" << highGCPls[idx][i][0][j]
                           << "\t     " << highGCPls[idx][i][1][j] << "\t     "
                           << highGCPls[idx][i][2][j] << endl;
            }

                ofstrm << "#\n# Medium Probabilities file for the exon model\n#\n"
                   << endl;


            ofstrm << "\n# Medium Die P_l's\n[Medium P_ls]" << endl;
            ofstrm << "# k = " << k << endl;
            for( int i = 0; i <= k; i++ ){
                ofstrm << "# l=\n" <<  i << endl;
                Seq2Int s2i(i+1);
                ofstrm << "# Values" << endl;
                for( int j = 0; j < mediumGCPls[idx][i][0].size(); j++ )
                    ofstrm << s2i.inv(j) << "\t" << mediumGCPls[idx][i][0][j]
                           << "\t     " << mediumGCPls[idx][i][1][j] << "\t     "
                           << mediumGCPls[idx][i][2][j] << endl;
            }

                ofstrm << "#\n# Low Probabilities file for the exon model\n#\n"
                   << endl;


            ofstrm << "\n# Low Die P_l's\n[Low P_ls]" << endl;
            ofstrm << "# k = " << k << endl;
            for( int i = 0; i <= k; i++ ){
                ofstrm << "# l=\n" <<  i << endl;
                Seq2Int s2i(i+1);
                ofstrm << "# Values" << endl;
                for( int j = 0; j < lowGCPls[idx][i][0].size(); j++ )
                    ofstrm << s2i.inv(j) << "\t" << lowGCPls[idx][i][0][j]
                           << "\t     " << lowGCPls[idx][i][1][j] << "\t     "
                           << lowGCPls[idx][i][2][j] << endl;
            }


	    /************* frequencies of the (k+1)-patterns

	    Seq2Int s2i(k+1);
	    ofstrm << "\n# normal patterns" << endl;
	    for( int j = 0; j < patterncount[0].size(); j++ )
		ofstrm << s2i.inv(j) << "\t" << patterncount[0][j]
		       << "\t     " << patterncount[1][j] << "\t     "
		       << patterncount[2][j] << endl;
	    
	    ofstrm << "\n# initial patterns" << endl;
	    for( int j = 0; j < initpatterncount[0].size(); j++ )
		ofstrm << s2i.inv(j) << "\t" << initpatterncount[0][j]
		       << "\t     " << initpatterncount[1][j] << "\t     "
		       << initpatterncount[2][j] << endl;
	    
	    ofstrm << "\n# internal exon terminal patterns" << endl;
	    for( int j = 0; j < etpatterncount[0].size(); j++ )
		ofstrm << s2i.inv(j) << "\t" << etpatterncount[0][j]
		       << "\t     " << etpatterncount[1][j] << "\t     "
		       << etpatterncount[2][j] << endl;

		       //  ************/


	    ofstrm << "\n# translation initiation motif\n[TRANSINIT]" << endl;
	    GCtransInitMotif[idx].write(ofstrm);

	    ofstrm << "\n# dss upstream motif, reading frame 0(reverse)\n[ETMOTIF0]" << endl;
	    GCetMotif[idx][0]->write(ofstrm);
	    ofstrm << "\n# dss upstream motif, reading frame 1(reverse)\n[ETMOTIF1]" << endl;
	    GCetMotif[idx][1]->write(ofstrm);
	    ofstrm << "\n# dss upstream motif, reading frame 2(reverse)\n[ETMOTIF2]" << endl;
	    GCetMotif[idx][2]->write(ofstrm);
	    

	    // MULTIMODEL
	    //ofstrm << "\n[CONTENT MODELS]" << endl;
	    //baumwelch->printAllContentModels(ofstrm);

	    // ENDE MULTIMODEL
/*
	    if (verbosity) 
		cout << "Writing amino acid patterns to file " << fname << "." << endl;
	    aadep.printTrans(ofstrm);
*/
	    // for human readability only: emission probabilities 

	    int precision = LLDouble::getOutputPrecision();
	    /*LLDouble::setOutputPrecision(4);
	    ofstrm << "\n\n#\n# Emission probabilities\n#\n[EMISSION]\n";
	    ofstrm << "# Size of vector\n" << GCemiprobs[idx].probs[0].size() << endl;
	    ofstrm << "# k : order of the markov model\n" << k << endl;   
	    ofstrm << "# patpseudocount (pseudocount of sequence patterns)\n" << patpseudo << endl;
	    ofstrm << "# Probabilities\n# Format: pattern win0 win1 win2\n";
	    for( int i = 0; i < GCemiprobs[idx].probs[0].size(); i++ ){
		ofstrm << s2i_e.inv(i);
		for (int f=0; f<3; f++)
		    ofstrm << "\t" << GCemiprobs[idx].probs[f][i];
		ofstrm << endl;
	    }
	    LLDouble::setOutputPrecision(precision);*/
        
	    LLDouble::setOutputPrecision(4);
        ofstrm << "\n\n#\n# High Emission probabilities\n#\n[High EMISSION]\n";
	    ofstrm << "# Size of vector\n" << highGCemiprobs[idx].probs[0].size() << endl;
	    ofstrm << "# k : order of the markov model\n" << k << endl;   
	    ofstrm << "# patpseudocount (pseudocount of sequence patterns)\n" << patpseudo << endl;
	    ofstrm << "# Probabilities\n# Format: pattern win0 win1 win2\n";
	    for( int i = 0; i < highGCemiprobs[idx].probs[0].size(); i++ ){
		ofstrm << s2i_e.inv(i);
		for (int f=0; f<3; f++)
		    ofstrm << "\t" << highGCemiprobs[idx].probs[f][i];
		ofstrm << endl;
	    }
	    LLDouble::setOutputPrecision(precision);
        
	    LLDouble::setOutputPrecision(4);
	    ofstrm << "\n\n#\n# Medium Emission probabilities\n#\n[Medium EMISSION]\n";
	    ofstrm << "# Size of vector\n" << mediumGCemiprobs[idx].probs[0].size() << endl;
	    ofstrm << "# k : order of the markov model\n" << k << endl;   
	    ofstrm << "# patpseudocount (pseudocount of sequence patterns)\n" << patpseudo << endl;
	    ofstrm << "# Probabilities\n# Format: pattern win0 win1 win2\n";
	    for( int i = 0; i < mediumGCemiprobs[idx].probs[0].size(); i++ ){
		ofstrm << s2i_e.inv(i);
		for (int f=0; f<3; f++)
		    ofstrm << "\t" << mediumGCemiprobs[idx].probs[f][i];
		ofstrm << endl;
	    }
	    LLDouble::setOutputPrecision(precision);        
        
	    LLDouble::setOutputPrecision(4);
	    ofstrm << "\n\n#\n# Low Emission probabilities\n#\n[Low EMISSION]\n";
	    ofstrm << "# Size of vector\n" << lowGCemiprobs[idx].probs[0].size() << endl;
	    ofstrm << "# k : order of the markov model\n" << k << endl;   
	    ofstrm << "# patpseudocount (pseudocount of sequence patterns)\n" << patpseudo << endl;
	    ofstrm << "# Probabilities\n# Format: pattern win0 win1 win2\n";
	    for( int i = 0; i < lowGCemiprobs[idx].probs[0].size(); i++ ){
		ofstrm << s2i_e.inv(i);
		for (int f=0; f<3; f++)
		    ofstrm << "\t" << lowGCemiprobs[idx].probs[f][i];
		ofstrm << endl;
	    }
	    LLDouble::setOutputPrecision(precision);        
	    
	    /*ofstrm << "\n\n#\n# Initial emission probabilities\n#\n[INITEMISSION]\n";
	    ofstrm << "# Size of vector\n" << GCinitemiprobs[idx][0].size() << endl;
	    ofstrm << "# k : order of the markov model\n" << k << endl;   
	    ofstrm << "# patpseudocount (pseudocount of sequence patterns)\n" << patpseudo << endl;
	    ofstrm << "# Probabilities\n# Format: pattern win0 win1 win2\n";
	    for( int i = 0; i < GCinitemiprobs[idx][0].size(); i++ ){
		ofstrm << s2i_e.inv(i);
		for (int f=0; f<3; f++)
		    ofstrm << "     \t" << GCinitemiprobs[idx][f][i];
		ofstrm << endl;
	    }*/
	    
	    ofstrm << "\n\n#\n# High Initial emission probabilities\n#\n[High INITEMISSION]\n";
	    ofstrm << "# Size of vector\n" << highGCinitemiprobs[idx][0].size() << endl;
	    ofstrm << "# k : order of the markov model\n" << k << endl;
	    ofstrm << "# patpseudocount (pseudocount of sequence patterns)\n" << patpseudo << endl;
	    ofstrm << "# Probabilities\n# Format: pattern win0 win1 win2\n";
	    for( int i = 0; i < highGCinitemiprobs[idx][0].size(); i++ ){
		ofstrm << s2i_e.inv(i);
		for (int f=0; f<3; f++)
		    ofstrm << "     \t" << highGCinitemiprobs[idx][f][i];
		ofstrm << endl;
	    }

	    ofstrm << "\n\n#\n# Medium Initial emission probabilities\n#\n[Medium INITEMISSION]\n";
	    ofstrm << "# Size of vector\n" << mediumGCinitemiprobs[idx][0].size() << endl;
	    ofstrm << "# k : order of the markov model\n" << k << endl;
	    ofstrm << "# patpseudocount (pseudocount of sequence patterns)\n" << patpseudo << endl;
	    ofstrm << "# Probabilities\n# Format: pattern win0 win1 win2\n";
	    for( int i = 0; i < mediumGCinitemiprobs[idx][0].size(); i++ ){
		ofstrm << s2i_e.inv(i);
		for (int f=0; f<3; f++)
		    ofstrm << "     \t" << mediumGCinitemiprobs[idx][f][i];
		ofstrm << endl;
	    }

	    ofstrm << "\n\n#\n# Low Initial emission probabilities\n#\n[Low INITEMISSION]\n";
	    ofstrm << "# Size of vector\n" << lowGCinitemiprobs[idx][0].size() << endl;
	    ofstrm << "# k : order of the markov model\n" << k << endl;
	    ofstrm << "# patpseudocount (pseudocount of sequence patterns)\n" << patpseudo << endl;
	    ofstrm << "# Probabilities\n# Format: pattern win0 win1 win2\n";
	    for( int i = 0; i < lowGCinitemiprobs[idx][0].size(); i++ ){
		ofstrm << s2i_e.inv(i);
		for (int f=0; f<3; f++)
		    ofstrm << "     \t" << lowGCinitemiprobs[idx][f][i];
		ofstrm << endl;
	    }

	    /*ofstrm << "\n\n#\n# Internal exon terminal emission probabilities\n#\n[ETEMISSION]\n";
	    ofstrm << "# Size of vector\n" << GCetemiprobs[idx][0].size() << endl;
	    ofstrm << "# k : order of the markov model\n" << k << endl;   
	    ofstrm << "# patpseudocount (pseudocount of sequence patterns)\n" << patpseudo << endl;
	    ofstrm << "# Probabilities\n# Format: pattern win0 win1 win2\n";
	    for( int i = 0; i < GCetemiprobs[idx][0].size(); i++ ){
		ofstrm << s2i_e.inv(i);
		for (int f=0; f<3; f++)
		    ofstrm << "     \t" << GCetemiprobs[idx][f][i];
		ofstrm << endl;
	    }*/

	    ofstrm << "\n\n#\n# High Internal exon terminal emission probabilities\n#\n[High ETEMISSION]\n";
	    ofstrm << "# Size of vector\n" << highGCetemiprobs[idx][0].size() << endl;
	    ofstrm << "# k : order of the markov model\n" << k << endl;
	    ofstrm << "# patpseudocount (pseudocount of sequence patterns)\n" << patpseudo << endl;
	    ofstrm << "# Probabilities\n# Format: pattern win0 win1 win2\n";
	    for( int i = 0; i < highGCetemiprobs[idx][0].size(); i++ ){
		ofstrm << s2i_e.inv(i);
		for (int f=0; f<3; f++)
		    ofstrm << "     \t" << highGCetemiprobs[idx][f][i];
		ofstrm << endl;
	    }

	    ofstrm << "\n\n#\n# Medium Internal exon terminal emission probabilities\n#\n[Medium ETEMISSION]\n";
	    ofstrm << "# Size of vector\n" << mediumGCetemiprobs[idx][0].size() << endl;
	    ofstrm << "# k : order of the markov model\n" << k << endl;
	    ofstrm << "# patpseudocount (pseudocount of sequence patterns)\n" << patpseudo << endl;
	    ofstrm << "# Probabilities\n# Format: pattern win0 win1 win2\n";
	    for( int i = 0; i < mediumGCetemiprobs[idx][0].size(); i++ ){
		ofstrm << s2i_e.inv(i);
		for (int f=0; f<3; f++)
		    ofstrm << "     \t" << mediumGCetemiprobs[idx][f][i];
		ofstrm << endl;
	    }

	    ofstrm << "\n\n#\n# Low Internal exon terminal emission probabilities\n#\n[Low ETEMISSION]\n";
	    ofstrm << "# Size of vector\n" << lowGCetemiprobs[idx][0].size() << endl;
	    ofstrm << "# k : order of the markov model\n" << k << endl;
	    ofstrm << "# patpseudocount (pseudocount of sequence patterns)\n" << patpseudo << endl;
	    ofstrm << "# Probabilities\n# Format: pattern win0 win1 win2\n";
	    for( int i = 0; i < lowGCetemiprobs[idx][0].size(); i++ ){
		ofstrm << s2i_e.inv(i);
		for (int f=0; f<3; f++)
		    ofstrm << "     \t" << lowGCetemiprobs[idx][f][i];
		ofstrm << endl;
	    }


	    /*
	    ofstrm << "# number of reading frames\n" << 3 << endl;
	    ofstrm << "# number of types\n" << numModels << endl;
	    ofstrm << "# probabilities of the emissions" << endl;
	    for (int i=0; i < emission[0].getRowSize(); i++) {
		ofstrm << s2i_e.inv(i) << "  ";
		for (int t=0; t < numModels; t++) {
		    for (int f=0; f < 3; f++) {
			ofstrm << setw(12) << emission[t][f][i];
		    }
		    if (t < numModels) {
			ofstrm << "  ";
		    }
		}
		ofstrm << endl;
	    }
	    */
	    ofstrm.close();
	}
	catch( ProjectError ) {
	    if (verbosity)    
		cout << "Exon model parameters not saved." << endl;
	}
}



/*
 * ===[ ExonModel::computeLengthDistributions]====================================
 */
void ExonModel::computeLengthDistributions( ){

    Smooth sm(minwindowcount, slope_of_bandwidth);
  
    /*
     * compute length distributions
     *
     */
    lenDistSingle.assign(max_exon_length+1, 0.0);
    lenDistInitial.assign(max_exon_length+1, 0.0);
    lenDistInternal.assign(max_exon_length+1, 0.0);
    lenDistTerminal.assign(max_exon_length+1, 0.0);
    
    highlenDistSingle.assign(max_exon_length+1, 0.0);
    highlenDistInitial.assign(max_exon_length+1, 0.0);
    highlenDistInternal.assign(max_exon_length+1, 0.0);
    highlenDistTerminal.assign(max_exon_length+1, 0.0);

    mediumlenDistSingle.assign(max_exon_length+1, 0.0);
    mediumlenDistInitial.assign(max_exon_length+1, 0.0);
    mediumlenDistInternal.assign(max_exon_length+1, 0.0);
    mediumlenDistTerminal.assign(max_exon_length+1, 0.0);

    lowlenDistSingle.assign(max_exon_length+1, 0.0);
    lowlenDistInitial.assign(max_exon_length+1, 0.0);
    lowlenDistInternal.assign(max_exon_length+1, 0.0);
    lowlenDistTerminal.assign(max_exon_length+1, 0.0);


    numHugeSingle   = numHugeExonsOfType[singleGHigh] + numHugeExonsOfType[singleGMedium] + numHugeExonsOfType[singleGLow];

    highnumHugeSingle   = numHugeExonsOfType[singleGHigh];
    mediumnumHugeSingle   = numHugeExonsOfType[singleGMedium];
    lownumHugeSingle   = numHugeExonsOfType[singleGLow];


    numHugeInitial  = numHugeExonsOfType[initial0High] + numHugeExonsOfType[initial0Medium] + numHugeExonsOfType[initial0Low]
                    +  numHugeExonsOfType[initial1High] + numHugeExonsOfType[initial1Medium] + numHugeExonsOfType[initial1Low]
                    + numHugeExonsOfType[initial2High] + numHugeExonsOfType[initial2Medium] + numHugeExonsOfType[initial2Low];

    highnumHugeInitial  = numHugeExonsOfType[initial0High]
                    +  numHugeExonsOfType[initial1High]
                    + numHugeExonsOfType[initial2High] ;
    mediumnumHugeInitial  = numHugeExonsOfType[initial0Medium]
                    + numHugeExonsOfType[initial1Medium]
                    + numHugeExonsOfType[initial2Medium] ;
    lownumHugeInitial  =  numHugeExonsOfType[initial0Low]
                    +  numHugeExonsOfType[initial1Low]
                    +  numHugeExonsOfType[initial2Low];

    numHugeInternal = numHugeExonsOfType[internal0High] + numHugeExonsOfType[internal0Medium] + numHugeExonsOfType[internal0Low]
                    +  numHugeExonsOfType[internal1High] + numHugeExonsOfType[internal1Medium] + numHugeExonsOfType[internal1Low]
	+ numHugeExonsOfType[internal2High] + numHugeExonsOfType[internal2Medium] + numHugeExonsOfType[internal2Low];

    highnumHugeInternal = numHugeExonsOfType[internal0High]
                        +  numHugeExonsOfType[internal1High]
    	+ numHugeExonsOfType[internal2High] ;

    mediumnumHugeInternal = numHugeExonsOfType[internal0Medium]
                        +  numHugeExonsOfType[internal1Medium]
    	+ numHugeExonsOfType[internal2Medium] ;

    lownumHugeInternal = numHugeExonsOfType[internal0Low]
                        +  numHugeExonsOfType[internal1Low]
    	+ numHugeExonsOfType[internal2Low] ;


    numHugeTerminal = numHugeExonsOfType[terminalHigh] + numHugeExonsOfType[terminalMedium] + numHugeExonsOfType[terminalLow];

    highnumHugeTerminal = numHugeExonsOfType[terminalHigh];
    mediumnumHugeTerminal = numHugeExonsOfType[terminalMedium] ;
    lownumHugeTerminal = numHugeExonsOfType[terminalLow];

    numSingle   = numExonsOfType[singleGHigh] + numExonsOfType[singleGMedium] + numExonsOfType[singleGLow];

    highnumSingle   = numExonsOfType[singleGHigh] ;
    mediumnumSingle   = numExonsOfType[singleGMedium];
    lownumSingle   = numExonsOfType[singleGLow];

    numInitial  = numExonsOfType[initial0High] + numExonsOfType[initial0Medium] + numExonsOfType[initial0Low]
                + numExonsOfType[initial1High] + numExonsOfType[initial1Medium] + numExonsOfType[initial1Low]
                + numExonsOfType[initial2High] + numExonsOfType[initial2Medium] + numExonsOfType[initial2Low];

    highnumInitial  = numExonsOfType[initial0High]
                    + numExonsOfType[initial1High]
                    + numExonsOfType[initial2High] ;

    mediumnumInitial  = numExonsOfType[initial0Medium]
                    + numExonsOfType[initial1Medium]
                    + numExonsOfType[initial2Medium] ;

    lownumInitial  = numExonsOfType[initial0Low]
                    + numExonsOfType[initial1Low]
                    + numExonsOfType[initial2Low] ;


    numInternal = numExonsOfType[internal0High] + numExonsOfType[internal0Medium] + numExonsOfType[internal0Low]
                + numExonsOfType[internal1High] + numExonsOfType[internal1Medium] + numExonsOfType[internal1Low]
                + numExonsOfType[internal2High] + numExonsOfType[internal2Medium] + numExonsOfType[internal2Low];

    highnumInternal = numExonsOfType[internal0High]
                    + numExonsOfType[internal1High]
                    + numExonsOfType[internal2High];

    mediumnumInternal = numExonsOfType[internal0Medium]
                    + numExonsOfType[internal1Medium]
                    + numExonsOfType[internal2Medium];

    lownumInternal = numExonsOfType[internal0Low]
                    + numExonsOfType[internal1Low]
                    + numExonsOfType[internal2Low];

    numTerminal = numExonsOfType[terminalHigh] + numExonsOfType[terminalMedium] + numExonsOfType[terminalLow];

    highnumTerminal = numExonsOfType[terminalHigh] ;
    mediumnumTerminal = numExonsOfType[terminalMedium] ;
    lownumTerminal = numExonsOfType[terminalLow];
   
    // compute mean exon lengths
    /*
    if (verbosity) {
	long int sum;
	int num, i;
	cout << "single, initial, internal, terminal mean exon lengths :" << endl;
	sum = 0; num = 0;
	for (i=0; i<= exonLenD; i++) {
	    sum += i*lenCountSingle[i]; num += lenCountSingle[i];
	}
	if (num>0)
	    cout << sum/num << "\t";
	else 
	    cout << "n.a." << "\t";
	sum = 0; num = 0;
	for (i=0; i<= exonLenD; i++) {
	    sum += i*lenCountInitial[i]; num += lenCountInitial[i];
	}
	if (num>0)
	    cout << sum/num << "\t";
	else 
	    cout << "n.a." << "\t";
	sum = 0; num = 0;
	for (i=0; i<= exonLenD; i++) {
	    sum += i*lenCountInternal[i]; num += lenCountInternal[i];
	}
	if (num>0)
	    cout << sum/num << "\t";
	else 
	    cout << "n.a." << "\t";
	sum = 0; num = 0;
	for (i=0; i<= exonLenD; i++) {
	    sum += i*lenCountTerminal[i]; num += lenCountTerminal[i];
	}
	if (num>0)
	    cout << sum/num << endl;
	else 
	    cout << "n.a." << endl;
    }
    */


/*
    if (numHugeSingle>0) {
	cerr << numHugeSingle << " single exons were "
	     << " longer than " << exonLenD << endl;
    }
    if (numHugeInitial>0) {
	cerr << numHugeInitial << " initial exons were "
	     << " longer than " << exonLenD << endl;
    }
    if (numHugeInternal>0) {
	cerr << numHugeInternal << " internal exons were "
	     << " longer than " << exonLenD << endl;
    } 
    if (numHugeTerminal>0) {
	cerr << numHugeTerminal << " terminal exons were "
	     << " longer than " << exonLenD << endl;
    }
*/
    sm.smoothCounts(lenCountSingle, lenDistSingle);

    sm.smoothCounts(highlenCountSingle, highlenDistSingle);
    sm.smoothCounts(mediumlenCountSingle, mediumlenDistSingle);
    sm.smoothCounts(lowlenCountSingle, lowlenDistSingle);

    sm.smoothCounts(lenCountInitial, lenDistInitial);

    sm.smoothCounts(highlenCountInitial, highlenDistInitial);
    sm.smoothCounts(mediumlenCountInitial, mediumlenDistInitial);
    sm.smoothCounts(lowlenCountInitial, lowlenDistInitial);

    sm.smoothCounts(lenCountInternal, lenDistInternal);

    sm.smoothCounts(highlenCountInternal, highlenDistInternal);
    sm.smoothCounts(mediumlenCountInternal, mediumlenDistInternal);
    sm.smoothCounts(lowlenCountInternal, lowlenDistInternal);

    sm.smoothCounts(lenCountTerminal, lenDistTerminal);

    sm.smoothCounts(highlenCountTerminal, highlenDistTerminal);
    sm.smoothCounts(mediumlenCountTerminal, mediumlenDistTerminal);
    sm.smoothCounts(lowlenCountTerminal, lowlenDistTerminal);

    // make single exons shorter than the minimal coding length impossible in the first place
    for (int i=0; i < Constant::min_coding_len && i < lenDistSingle.size(); i++)
	lenDistSingle[i] = 0.0;

    for (int i=0; i < Constant::min_coding_len && i < highlenDistSingle.size(); i++)
	highlenDistSingle[i] = 0.0;

    for (int i=0; i < Constant::min_coding_len && i < mediumlenDistSingle.size(); i++)
	mediumlenDistSingle[i] = 0.0;

    for (int i=0; i < Constant::min_coding_len && i < lowlenDistSingle.size(); i++)
	lowlenDistSingle[i] = 0.0;

    /*
    scaleDblVector(lenDistSingle, Double(numSingle-numHugeSingle) / numSingle);

    scaleDblVector(highlenDistSingle, Double(numSingle-numHugeSingle) / numSingle);
    scaleDblVector(mediumlenDistSingle, Double(numSingle-numHugeSingle) / numSingle);
    scaleDblVector(lowlenDistSingle, Double(numSingle-numHugeSingle) / numSingle);

    scaleDblVector(lenDistInitial, Double(numInitial-numHugeInitial) / numInitial);

    scaleDblVector(highlenDistInitial, Double(numInitial-numHugeInitial) / numInitial);
    scaleDblVector(mediumlenDistInitial, Double(numInitial-numHugeInitial) / numInitial);
    scaleDblVector(lowlenDistInitial, Double(numInitial-numHugeInitial) / numInitial);

    scaleDblVector(lenDistInternal, Double(numInternal-numHugeInternal) / numInternal);
    scaleDblVector(lenDistTerminal, Double(numTerminal-numHugeTerminal) / numTerminal);

    scaleDblVector(highlenDistTerminal, Double(numTerminal-numHugeTerminal) / numTerminal);
    scaleDblVector(mediumlenDistTerminal, Double(numTerminal-numHugeTerminal) / numTerminal);
    scaleDblVector(lowlenDistTerminal, Double(numTerminal-numHugeTerminal) / numTerminal);
    */

    scaleDblVector(lenDistSingle, Double(numSingle-numHugeSingle) / numSingle);

	scaleDblVector(highlenDistSingle, Double(highnumSingle-highnumHugeSingle) / highnumSingle);
	scaleDblVector(mediumlenDistSingle, Double(mediumnumSingle-mediumnumHugeSingle) / mediumnumSingle);
	scaleDblVector(lowlenDistSingle, Double(lownumSingle-lownumHugeSingle) / lownumSingle);

	scaleDblVector(lenDistInitial, Double(numInitial-numHugeInitial) / numInitial);

	scaleDblVector(highlenDistInitial, Double(highnumInitial-highnumHugeInitial) / highnumInitial);
	scaleDblVector(mediumlenDistInitial, Double(mediumnumInitial-mediumnumHugeInitial) / mediumnumInitial);
	scaleDblVector(lowlenDistInitial, Double(lownumInitial-lownumHugeInitial) / lownumInitial);

	scaleDblVector(lenDistInternal, Double(numInternal-numHugeInternal) / numInternal);

	scaleDblVector(highlenDistInternal, Double(highnumInternal-highnumHugeInternal) / highnumInternal);
	scaleDblVector(mediumlenDistInternal, Double(mediumnumInternal-mediumnumHugeInternal) / mediumnumInternal);
	scaleDblVector(lowlenDistInternal, Double(lownumInternal-lownumHugeInternal) / lownumInternal);

	scaleDblVector(lenDistTerminal, Double(numTerminal-numHugeTerminal) / numTerminal);

	scaleDblVector(highlenDistTerminal, Double(highnumTerminal-highnumHugeTerminal) / highnumTerminal);
	scaleDblVector(mediumlenDistTerminal, Double(mediumnumTerminal-mediumnumHugeTerminal) / mediumnumTerminal);
	scaleDblVector(lowlenDistTerminal, Double(lownumTerminal-lownumHugeTerminal) / lownumTerminal);


    fillTailsOfLengthDistributions();
    fillTailsOfLengthDistributionsHigh();
    fillTailsOfLengthDistributionsMedium();
    fillTailsOfLengthDistributionsLow();
}



/*
 * ===[ ExonModel::processExons ]=========================================
 */
void ExonModel::processExons( const Gene* gene ){
//	cout << "sequence: " << gene->seqname << endl;
    State   *exon = gene->exons;
    curwin = 0;
    if( !exon->next ){
	try{
	    processSingleExon( exon );
	} catch (ExonModelError e) {
	    if (verbosity)
	      cerr << "gene " << gene->geneid << " transcr. " << gene->id << " in sequence " << gene->seqname << ": " << e.getMessage() << endl;
	}
    }
    else{
	try{
	    processInitialExon( exon );
	} catch( ExonModelError e ){
	    if (verbosity)
	      cerr << "gene " << gene->geneid << " transcr. " << gene->id << " in sequence " << gene->seqname << ": " << e.getMessage() << endl;
	}
	exon = exon->next;
	//prewin=-1;
	while( exon->next ){
	    try{
		processInternalExon( exon );
	    } catch( ExonModelError e ){
	      if (verbosity)
		    cerr << "gene " << gene->geneid << " transcr. " << gene->id << " in sequence " << gene->seqname << ": " << e.getMessage() << endl;
	    }
	    exon = exon->next;
	} 
	try{
	    processTerminalExon( exon );
	} catch( ExonModelError e ){
	    if (verbosity)
		cerr << "gene " << gene->geneid << " transcr. " << gene->id << " in sequence " << gene->seqname << ": " << e.getMessage() << endl;
	}
    }
}


/*
 * ===[ ExonModel::processSingleExon ]====================================
 */
void ExonModel::processSingleExon( State* exon) {
    int length = exon->end-exon->begin+1; 
    if (length < STARTCODON_LEN + STOPCODON_LEN) {
	throw ExonModelError("Single training exon too short.");
    } 
    if (!onStart(sequence+exon->begin)) {
      string msg("Single exon gene doesn't begin with atg codon but with ");
      msg += string(sequence+exon->begin, 3);
	throw ExonModelError(msg);
    }
   
    /*
    char motif[trans_init_window+1];
    strncpy(motif, sequence + exon->begin-12, trans_init_window);
    motif[trans_init_window] = '\0';
    cout << "trans init single " << motif << endl;
    */
    if (length > trans_init_window && exon->begin >= trans_init_window + k)
	transInitMotif->addSequence(sequence + exon->begin - trans_init_window, gweight);

    curwin = STARTCODON_LEN + k;

    const char *beginOfInnerSeq, *endOfInnerSeq, *endOfInit;
    beginOfInnerSeq =  sequence + exon->begin + STARTCODON_LEN + k;
    endOfInnerSeq = sequence + exon->end - STOPCODON_LEN;
    endOfInit = beginOfInnerSeq + Constant::init_coding_len - 1;

    processInnerSequence(beginOfInnerSeq, endOfInit, 2); // init part
    processInnerSequence(endOfInit+1, endOfInnerSeq, 0); // normal 

    int stppos = exon->end - STOPCODON_LEN + 1;
    if (strncmp(sequence + stppos, "taa", 3)==0)
	ochrecount++;
    else if (strncmp(sequence + stppos, "tag", 3)==0)
	ambercount++;
    else if (strncmp(sequence + stppos, "tga", 3)==0)
	opalcount++;
    else 
	throw ExonModelError("Single exon doesn't end in stop codon. Variable stopCodonExcludedFromCDS set right?");
    
    if (!hasLenDist) {
    	int na=0, nc=0, ng=0, nt=0;
    	    for (const char* cur = sequence+exon->begin; cur <= sequence+exon->end; cur++) {
    	        switch (*cur) {
    	            case 'a':
    	                na++;
    	                break;
    	            case 'c':
    	                nc++;
    	                break;
    	            case 'g':
    	                ng++;
    	                break;
    	            case 't':
    	                nt++;
    	                break;
    	        }
    	    }
    	    int ngc = ng + nc;
    	    double rgc = ngc / ((na+nc+ng+nt)*1.0);
    	    int type = -1;
    	    if (rgc < Properties::getdoubleProperty("lowT"))
    	        type = 3;
    	    else if (rgc >= Properties::getdoubleProperty("lowT") && rgc < Properties::getdoubleProperty("highT"))
    	        type = 2;
    	    else if (rgc >= Properties::getdoubleProperty("highT"))
    	        type = 1;
    	    if (exon->type == singleGHigh) {
    	    	switch(type){
    	    		case 1:
    					exon->type = singleGHigh;
    					break;
    				case 2:
    					exon->type = singleGMedium;
    					break;
    				case 3:
    					exon->type = singleGLow;
    					break;
    				default:
    					break;
    	    	}
    	    } else if (exon->type == rsingleGHigh){
    	    	switch(type){
    	    		case 1:
    					exon->type = rsingleGHigh;
    					break;
    				case 2:
    					exon->type = rsingleGMedium;
    					break;
    				case 3:
    					exon->type = rsingleGLow;
    					break;
    				default:
    					break;
    	    	}
    	    }
            cout << "single exon: "; exon->print();
            cout << "type " << type << endl;
            cout << "rgc: " << rgc << endl;
      if (length <= exonLenD) {
	lenCountSingle[length] += 1; // gweight;

	    if (type==3) {
	    	lowlenCountSingle[length] += 1; // gweight
	    } else if (type==2) {
	    	mediumlenCountSingle[length] += 1; // gweight;
	    } else if (type==1) {
	    	highlenCountSingle[length] += 1; // gweight;
	    }
      } else {
	numHugeExonsOfType[singleGHigh+type-1] += 1; // gweight;
      }
      numExonsOfType[singleGHigh+type-1] += 1; // gweight ;
    }
}


/*
 * ===[ ExonModel::processInitialExon ]=====================================
 */
void ExonModel::processInitialExon( State* exon){
    int oldwin = curwin;
    int length = exon->end-exon->begin+1; 
    

    if( exon->end - exon->begin + 1 < STARTCODON_LEN )
    {
	curwin = curwin + length;
	throw ExonModelError( "Initial exon has length < 3!" );
    }
    if (!onStart(sequence+exon->begin)) {
	curwin = curwin + length;
	throw ExonModelError("Initial Exon doesn't begin with start codon.");
    } 

    if (length > trans_init_window && exon->begin >= trans_init_window + k)
	transInitMotif->addSequence(sequence + exon->begin - trans_init_window, gweight);
    
    if (exon->end - Constant::dss_start - Constant::et_coding_len + 1 >=0 )
	etMotif[mod3(oldwin + length - Constant::dss_start - Constant::et_coding_len)]->addSequence(
	    sequence + exon->end - Constant::dss_start - Constant::et_coding_len + 1, gweight, true);

    /* test output for seqlogo/pictogram
    static char *fenster = new char[Constant::et_coding_len+1];
    strncpy( fenster, sequence + exon->end - Constant::dss_start - Constant::et_coding_len + 1, Constant::et_coding_len+1 );
    fenster[Constant::et_coding_len]='\0';
    cout << "initial  etMotif " << mod3(oldwin + length - Constant::dss_start - Constant::et_coding_len) << " " <<fenster << endl;
    */


    curwin = curwin + STARTCODON_LEN + k;
    try {
	const char *beginOfInnerSeq, *endOfInit, *endOfInnerSeq;
	beginOfInnerSeq = sequence + exon->begin + STARTCODON_LEN + k;
	endOfInnerSeq = sequence + exon->end - Constant::dss_start;
	endOfInit = beginOfInnerSeq + Constant::init_coding_len - 1;
	if (endOfInit > endOfInnerSeq) // at most up to the end of the sequence
	    endOfInit = endOfInnerSeq;
	
	processInnerSequence(beginOfInnerSeq, endOfInit, 2); // init part
	processInnerSequence(endOfInit+1, endOfInnerSeq);    // normal

    } catch( ExonModelError& e ){
	cerr << "ExonModel::processInitialExon: " << e.getMessage() << endl;
	throw e;
    }
    curwin = oldwin + length;

    int na=0, nc=0, ng=0, nt=0;
    for (const char* cur = sequence+exon->begin; cur <= sequence+exon->end; cur++) {
        switch (*cur) {
            case 'a':
                na++;
                break;
            case 'c':
                nc++;
                break;
            case 'g':
                ng++;
                break;
            case 't':
                nt++;
                break;
        }
    }
    int ngc = ng + nc;
    double rgc = ngc / ((na+nc+ng+nt)*1.0);
    int type = -1;
    if (rgc < Properties::getdoubleProperty("lowT"))
        type = 3;
    else if (rgc >= Properties::getdoubleProperty("lowT") && rgc < Properties::getdoubleProperty("highT"))
        type = 2;
    else if (rgc >= Properties::getdoubleProperty("highT"))
        type = 1;
    if (exon->type == initial0High) {
    	switch(type){
    		case 1:
				exon->type = initial0High;
				break;
			case 2:
				exon->type = initial0Medium;
				break;
			case 3:
				exon->type = initial0Low;
				break;
			default:
				break;
    	}
    } else if (exon->type == initial1High) {
    	switch(type){
    		case 1:
				exon->type = initial1High;
				break;
			case 2:
				exon->type = initial1Medium;
				break;
			case 3:
				exon->type = initial1Low;
				break;
			default:
				break;
    	}
    }else if (exon->type == initial2High) {
    	switch(type){
    		case 1:
				exon->type = initial2High;
				break;
			case 2:
				exon->type = initial2Medium;
				break;
			case 3:
				exon->type = initial2Low;
				break;
			default:
				break;
    	}
    }else if (exon->type == rinitialHigh) {
    	switch(type){
    		case 1:
				exon->type = rinitialHigh;
				break;
			case 2:
				exon->type = rinitialMedium;
				break;
			case 3:
				exon->type = rinitialLow;
				break;
			default:
				break;
    	}
    }
//    cout << "initial exon: "; exon->print();
//    cout << "type " << type << endl;
//    cout << "rgc: " << rgc << endl;

    if (!hasLenDist) {
	if (length <= exonLenD) {
	    lenCountInitial[length] += 1; // gweight;

	    if (type==3) {
	    	lowlenCountInitial[length] += 1;
	    } else if (type==2) {
	    	mediumlenCountInitial[length] += 1; // gweight;
	    } else if (type==1) {
	    	highlenCountInitial[length] += 1; // gweight;
	    }

	} else {
		if (type==1) {
			numHugeExonsOfType[initialExonHigh(mod3(curwin))] += 1; //gweight;
		} else if (type==2) {
			numHugeExonsOfType[initialExonMedium(mod3(curwin))] += 1; //gweight;
		} else if (type==3) {
			numHugeExonsOfType[initialExonLow(mod3(curwin))] += 1; //gweight;
		}
	}  
		if (type==1) {
			numExonsOfType[initialExonHigh(mod3(curwin))] += 1; //gweight;
		} else if (type==2) {
			numExonsOfType[initialExonMedium(mod3(curwin))] += 1; //gweight;
		} else if (type==3) {
			numExonsOfType[initialExonLow(mod3(curwin))] += 1; //gweight;
		}
    }
}


/*
 * ===[ ExonModel::processInternalExon ]=====================================
 *
 * cumcodinglength is the coding length of all exons before this one and is needed 
 * for training the init content models
 */
void ExonModel::processInternalExon( State* exon){
    int oldwin = curwin;
    int length = exon->end-exon->begin+1;
    
    curwin += Constant::ass_end + k;  
    try{

	const char *beginOfInnerSeq, *beginOfET, *endOfInnerSeq;
	beginOfInnerSeq = sequence + exon->begin + Constant::ass_end + k;
	endOfInnerSeq = sequence + exon->end - Constant::dss_start;
	beginOfET = endOfInnerSeq - Constant::et_coding_len + 1;
	if (beginOfET < beginOfInnerSeq)
	    beginOfET = beginOfInnerSeq;
	
	processInnerSequence(beginOfInnerSeq, beginOfET - 1);  // normal
	processInnerSequence(beginOfET, endOfInnerSeq, 3);     // internal exon term part

	etMotif[mod3(oldwin + length - Constant::dss_start - Constant::et_coding_len)]->addSequence(
	    sequence + exon->end - Constant::dss_start - Constant::et_coding_len + 1, gweight, true);
	/* test output for seqlogo/pictogram
	static char *fenster = new char[Constant::et_coding_len+1];
	strncpy( fenster, sequence + exon->end - Constant::dss_start - Constant::et_coding_len + 1, Constant::et_coding_len+1 );
	fenster[Constant::et_coding_len]='\0';
	cout << "internal etMotif " << mod3(oldwin + length - Constant::dss_start - Constant::et_coding_len) << " " <<fenster << endl;
	*/
	
    } catch( ExonModelError e ){
	cerr << "ExonModel::processInternalExon: " << e.getMessage() << endl;
	throw e;
    }
    
    curwin = oldwin + length;
    //transition of phases: cout << "internal " << mod3(oldwin) << "->" << mod3(curwin) << endl;
    
/* for determining transition probabilities
    if (prewin>=0) {
	cout << prewin << curwin%3 << endl;
    } 
    prewin=curwin%3;*/

    int na=0, nc=0, ng=0, nt=0;
            for (const char* cur = sequence+exon->begin; cur <= sequence+exon->end; cur++) {
                switch (*cur) {
                    case 'a':
                        na++;
                        break;
                    case 'c':
                        nc++;
                        break;
                    case 'g':
                        ng++;
                        break;
                    case 't':
                        nt++;
                        break;
                }
            }
            int ngc = ng + nc;
            double rgc = ngc / ((na+nc+ng+nt)*1.0);
            int type = -1;
    	    if (rgc < Properties::getdoubleProperty("lowT"))
    	        type = 3;
    	    else if (rgc >= Properties::getdoubleProperty("lowT") && rgc < Properties::getdoubleProperty("highT"))
    	        type = 2;
    	    else if (rgc >= Properties::getdoubleProperty("highT"))
    	        type = 1;
    	    if (exon->type == internal0High) {
    	    	switch(type){
    	    		case 1:
    					exon->type = internal0High;
    					break;
    				case 2:
    					exon->type = internal0Medium;
    					break;
    				case 3:
    					exon->type = internal0Low;
    					break;
    				default:
    					break;
    	    	}
    	    } else if (exon->type == internal1High){
    	    	switch(type){
    	    		case 1:
    					exon->type = internal1High;
    					break;
    				case 2:
    					exon->type = internal1Medium;
    					break;
    				case 3:
    					exon->type = internal1Low;
    					break;
    				default:
    					break;
    	    	}
    	    }else if (exon->type == internal2High) {
    				switch(type){
    					case 1:
								exon->type = internal2High;
								break;
							case 2:
								exon->type = internal2Medium;
								break;
							case 3:
								exon->type = internal2Low;
								break;
							default:
								break;
    				}
    			}else if (exon->type == rinternal0High) {
    				switch(type){
    					case 1:
								exon->type = rinternal0High;
								break;
							case 2:
								exon->type = rinternal0Medium;
								break;
							case 3:
								exon->type = rinternal0Low;
								break;
							default:
								break;
    				}
    			}else if (exon->type == rinternal1High) {
    				switch(type){
    					case 1:
								exon->type = rinternal1High;
								break;
							case 2:
								exon->type = rinternal1Medium;
								break;
							case 3:
								exon->type = rinternal1Low;
								break;
							default:
								break;
    				}
    			}else if (exon->type == rinternal2High) {
    				switch(type){
    					case 1:
								exon->type = rinternal2High;
								break;
							case 2:
								exon->type = rinternal2Medium;
								break;
							case 3:
								exon->type = rinternal2Low;
								break;
							default:
								break;
    				}
    			}
    			
//            cout << "internal exon: "; exon->print();
//            cout << "type " << type << endl;
//            cout << "rgc: " << rgc << endl;

    if (!hasLenDist) {
	if (length <= exonLenD) {
	    lenCountInternal[length] += 1; // gweight;
	    if (type==1) {
	    	highlenCountInternal[length] += 1; // gweight;
	    } else if (type==2) {
	    	mediumlenCountInternal[length] += 1; // gweight;
	    } else if (type==3) {
	    	lowlenCountInternal[length] += 1; // gweight;
	    }

	} else {
		if (type==1) {
			numHugeExonsOfType[internalExonHigh(mod3(curwin))] += 1; //gweight;
		} else if (type==2) {
			numHugeExonsOfType[internalExonMedium(mod3(curwin))] += 1; //gweight;
		} else if (type==3) {
			numHugeExonsOfType[internalExonLow(mod3(curwin))] += 1; //gweight;
		}
	}
	if (type==1) {
			numExonsOfType[internalExonHigh(mod3(curwin))] += 1; //gweight;
		}else if (type==2) {
			numExonsOfType[internalExonMedium(mod3(curwin))] += 1; //gweight;
		} else if (type==3) {
			numExonsOfType[internalExonLow(mod3(curwin))] += 1; //gweight;
		}
    }
}


/*
 * ===[ ExonModel::processTerminalExon ]=======================================
 */
void ExonModel::processTerminalExon( State* exon){
    int length = exon->end-exon->begin+1; 
    const char *beginOfInnerSeq, *endOfInnerSeq;

    curwin += Constant::ass_end + k;
    
    beginOfInnerSeq = sequence + exon->begin + Constant::ass_end + k;
    endOfInnerSeq = sequence + exon->end - STOPCODON_LEN;

    processInnerSequence( beginOfInnerSeq, endOfInnerSeq);

    int stppos = exon->end - STOPCODON_LEN + 1;
    if (strncmp(sequence + stppos, "taa", 3)==0)
	ochrecount++;
    else if (strncmp(sequence + stppos, "tag", 3)==0)
	ambercount++;
    else if (strncmp(sequence + stppos, "tga", 3)==0)
	opalcount++;
    else 
	throw ExonModelError("Terminal exon doesn't end in stop codon. Variable stopCodonExcludedFromCDS set right?");

    int na=0, nc=0, ng=0, nt=0;
        for (const char* cur = sequence+exon->begin; cur <= sequence+exon->end; cur++) {
            switch (*cur) {
                case 'a':
                    na++;
                    break;
                case 'c':
                    nc++;
                    break;
                case 'g':
                    ng++;
                    break;
                case 't':
                    nt++;
                    break;
            }
        }
        int ngc = ng + nc;
        double rgc = ngc / ((na+nc+ng+nt)*1.0);
        int type = -1;
	    if (rgc < Properties::getdoubleProperty("lowT"))
	        type = 3;
	    else if (rgc >= Properties::getdoubleProperty("lowT") && rgc < Properties::getdoubleProperty("highT"))
	        type = 2;
	    else if (rgc >= Properties::getdoubleProperty("highT"))
	        type = 1;
	    if (exon->type == terminalHigh) {
	    	switch(type){
	    		case 1:
					exon->type = terminalHigh;
					break;
				case 2:
					exon->type = terminalMedium;
					break;
				case 3:
					exon->type = terminalLow;
					break;
				default:
					break;
	    	}
	    } else if (exon->type == rterminal0High) {
	    	switch(type){
	    		case 1:
						exon->type = rterminal0High;
						break;
					case 2:
						exon->type = rterminal0Medium;
						break;
					case 3:
						exon->type = rterminal0Low;
						break;
					default:
						break;
	    	}
	    }else if (exon->type == rterminal1High) {
    				switch(type){
    					case 1:
								exon->type = rterminal1High;
								break;
							case 2:
								exon->type = rterminal1Medium;
								break;
							case 3:
								exon->type = rterminal1Low;
								break;
							default:
								break;
    				}
    			}else if (exon->type == rterminal2High) {
    				switch(type){
    					case 1:
								exon->type = rterminal2High;
								break;
							case 2:
								exon->type = rterminal2Medium;
								break;
							case 3:
								exon->type = rterminal2Low;
								break;
							default:
								break;
    				}
    			}
//        cout << "terminal exon: "; exon->print();
//        cout << "type " << type << endl;
//        cout << "rgc: " << rgc << endl;

    if (!hasLenDist) {
	if (length <= exonLenD) {
	    lenCountTerminal[length] += 1; // gweight;

	    if (type==1) {
	    	highlenCountTerminal[length] += 1; // gweight;
	    } else if (type==2) {
	    	mediumlenCountTerminal[length] += 1; // gweight;
	    } else if (type==3) {
	    	lowlenCountTerminal[length] += 1; // gweight;
	    }
	} else {
	    numHugeExonsOfType[terminalHigh+type-1] += 1; // gweight;
	}
	numExonsOfType[terminalHigh+type-1] += 1; // gweight;
    }
}


/*
 * ===[ ExonModel::processInnerSequence ]=================================
 * 'begin' is the first nucleotide, of which the emission is counted
 * 'end' is the last
 */
void ExonModel::processInnerSequence(const char *begin, const char* end, int modeltype){
    if (begin > end)
	return;
    Seq2Int s2i(k+1);
    
    int na=0, nc=0, ng=0, nt=0;
    for (const char* cur = begin; cur <= end; cur++) {
        switch (*cur) {
            case 'a':
                na++;
            case 'c':
                nc++;
            case 'g':
                ng++;
            case 't':
                nt++;
        }
    }
    int ngc = ng + nc;    
    double rgc = ngc / ((na+nc+ng+nt)*1.0);
    int type = -1;        
    if (rgc < Properties::getdoubleProperty("lowT"))
        type = 3;
    else if (rgc >= Properties::getdoubleProperty("lowT") && rgc < Properties::getdoubleProperty("highT"))
        type = 2;
    else if (rgc >= Properties::getdoubleProperty("highT"))
        type = 1;
    
    for( ; begin <= end; begin++){        
	curwin %= 3;
	if (curwin==0 && end-begin>=STOPCODON_LEN && GeneticCode::isStopcodon(begin)){
	    throw ExonModelError("in-frame stop codon");
	}
	try {
	    int pn = s2i(begin-k);
	    if (modeltype == 0){
            patterncount[ curwin ][pn] += gweight;           
            
            switch (type) {
                case 1:
                    patternhighcount[ curwin ][pn] += gweight;
                    break;
                case 2:
                    patternmediumcount[ curwin ][pn] += gweight;
                    break;
                case 3:
                    patternlowcount[ curwin ][pn] += gweight;
                    break;
                default:
                    break;
            }
        }
	    else if (modeltype == 2 ) {// initial part of the coding region
		initpatterncount[ curwin ][pn] += gweight;
			switch (type) {
						case 1:
							highinitpatterncount[ curwin ][pn] += gweight;
							break;
						case 2:
							mediuminitpatterncount[ curwin ][pn] += gweight;
							break;
						case 3:
							lowinitpatterncount[ curwin ][pn] += gweight;
							break;
						default:
							break;
			}
	    }
	    else if (modeltype == 3 ) // terminal part of an internal exon
		{
	    	etpatterncount[ curwin ][pn] += gweight;
	    	switch (type) {
	    	case 1:
	    		highetpatterncount[ curwin ][pn] += gweight;
	    		break;
	    	case 2:
	    		mediumetpatterncount[ curwin ][pn] += gweight;
	    		break;
	    	case 3:
	    		lowetpatterncount[ curwin ][pn] += gweight;
	    		break;
	    	default:
	    		break;
	    	}
		}
	    else 
		throw ProjectError("ExonModel::processInnerSequence: invalid model type");
	    gesbasen[curwin] += gweight;
	} catch (InvalidNucleotideError e) {}
	curwin++;
    }
}



/*
 * ===[ ExonModel::computeCMFromBW ]===================================
 *

void ExonModel::computeCMFromBW(){
     *
     * Emission Probabilities for order k emissions only
     * They are computed using the content models from the baumwelch object
     *
    Matrix< vector<Double> > emission(numModels, 3);

    for (int t=0; t<numModels; t++) {
    	ContentModel cm = baumwelch->getContentModel(t);
	for (int f=0; f<3; f++) {
	    emission[t][f].resize(cm.getNumPatterns());
	    computeEmiFromPat(cm.copyPatProb(f), emission[t][f], k);
	}
    }

     *
     * Initialize modelStartProbs - the vector of probabilities of the
     * different content models.
     *
    if (modelStartProbs)
	delete [] modelStartProbs;
    modelStartProbs = new Double[numModels];
    for (int t=0; t<numModels; t++) {
	modelStartProbs[t] = baumwelch->getModelTypeProb(t);
    }
}
*/
