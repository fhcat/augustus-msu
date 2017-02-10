/*****************************************************************************\
 * Filename : exonmodel.hh
 * Author   : Mario Stanke
 * Project  : HMM
 *
 *
 * Description: Training the exon model and implementing the algorithms.
 *
 *
 * Date       |   Author              |  Changes
 *------------|-----------------------|------------------------------------------
 * ?          | Stafilarakis          | creations of the class
 * 06.05.2002 | Stanke                | debugging and testing
 * 18.05.2002 | Stanke                | added length distributions
 * 24.03.2003 | Mario Stanke          | introducing shadow states
 * 31.08.2005 | Mario Stanke          | sequence is no longer a member of gene
\******************************************************************************/

#ifndef _EXONMODEL_HH
#define _EXONMODEL_HH

#include "statemodel.hh"


/*
 * The reading frame of an exon is the position of the nucleotide following the exon
 * in its codon starting counting at 0:
 * forward strand:
 *  *** | ***        index 0
 *  *** * | ** ***   index 1
 *  *** ** | * ***   index 2
 * reverse strand:
 *  *** | ***        index 2
 *  *** * | ** ***   index 1
 *  *** ** | * ***   index 0
 */

class ExonModelError : public ProjectError {
public:
    ExonModelError(string msg) : ProjectError(msg) {}
};


/*
 * Open Reading Frame (ORF), part of the exon model
 */
class OpenReadingFrame{
public:
    OpenReadingFrame(const char* dna, int _max_exon_length, int _n);
    OpenReadingFrame() {}
    ~OpenReadingFrame() {}
    int leftmostExonBegin(int frame, int base, bool forward);

private:
    // used internally, these call GeneticCode
    static bool isStopcodon(const char* dna);
    static bool isRCStopcodon(const char* dna);

    vector<Integer> nearestStopForward;
    vector<Integer> nearestStopReverse;
    int n;
    int max_exon_length;
    // static bool amber, ochre, opal; // now taken care of by GeneticCode
};

/*
 * ExonModel
 */
class ExonModel : public StateModel{
public:
    ExonModel();
    ~ExonModel();

    StateType getStateType() const {
	return etype;
    }
    void buildModel         ( const AnnoSequence* annoseq, int parIndex );
    void registerPars       ( Parameters* parameters);
    void printProbabilities ( int parIndex, BaseCount *bc, const char* suffix = NULL );
    void viterbiForwardAndSampling(ViterbiMatrixType&, ViterbiMatrixType&, int, int,
				   AlgorithmVariant, OptionListItem&) const;
    Double emiProbUnderModel(int begin, int end) const;
    Double endPartEmiProb(int end) const;
    Double notEndPartEmiProb(int beginOfStart, int right, int frameOfRight, Feature *exonparts) const;
    void initAlgorithms(Matrix<Double>&, int);
    static void storeGCPars(int idx);
	
    // class functions
    static void init();
    static void resetPars() {
	initAlgorithmsCalled = false;
    }
    static void readProbabilities(int parIndex);
    static void readAllParameters();
    static void resetModelCount(){exoncount = 0;};
    static int getMaxStateLen() { return max_exon_length + trans_init_window; }
    static void setORF() {
      if (orf)
	delete orf;
      orf = new OpenReadingFrame(sequence, max_exon_length, dnalen);
      initAlgorithmsCalled = false;
    }

private:
    void processExons(const Gene* gene);
    void processSingleExon(State* exon);
    void processInitialExon(State* exon);
    void processInternalExon(State* exon);
    void processTerminalExon(State* exon);
    void processInnerSequence(const char* begin, const char* end, int modeltype = 0);
    void buildProbabilities ( );
    Double seqProb(int endOfStart, int right, int frameOfRight) const;
    Double eTermSeqProb(int left, int right, int frameOfRight) const;
    Double initialSeqProb(int left, int right, int frameOfRight) const;
    void computeLowerOrderPats();
    void computeCMFromBW    ();
 
    // internal class functions
    static void computeLengthDistributions( );
    static void fillTailsOfLengthDistributions( );
    static void fillTailsOfLengthDistributionsHigh( );
    static void fillTailsOfLengthDistributionsMedium( );
    static void fillTailsOfLengthDistributionsLow( );
    static bool isHighState(StateType type);
    static bool isMediumState(StateType type);
    static bool isLowState(StateType type);

    StateType        etype;
    Integer          win,                 // reading frame of this state (fixed)
	             curwin;              // current reading frame during training
    int              beginPartLen;        // constant for every etype: endOfPred = beginOfStart - beginPartLen -1
    int              innerPartOffset;     //                           beginOfStart = beginOfBioExon + innerPartOffset
    int              innerPartEndOffset;  //                           right = endOfBioExon - innerPartEndOffset
    int              baseOffset;          //                           base = endOfBioExon - baseOffset
    Integer          gweight;
    // int prewin; 

    // class variables
    static vector<Integer> patterncount[3];     // {0,1,2}x{acgt}^(k+1), the reading frame is 
                                                // the position of the emitted (last) nucleotide
    static vector<Integer> patternhighcount[3];
    static vector<Integer> patternmediumcount[3];
    static vector<Integer> patternlowcount[3];

    static vector<Integer> initpatterncount[3]; // like above, just the first nucleotides of a gene

    static vector<Integer> highinitpatterncount[3]; // like above, just the first nucleotides of a gene
    static vector<Integer> mediuminitpatterncount[3]; // like above, just the first nucleotides of a gene
    static vector<Integer> lowinitpatterncount[3]; // like above, just the first nucleotides of a gene

    static vector<Integer> etpatterncount[3] ;  // internal exon ternimal part

    static vector<Integer> highetpatterncount[3] ;  // internal exon ternimal part
    static vector<Integer> mediumetpatterncount[3] ;  // internal exon ternimal part
    static vector<Integer> lowetpatterncount[3] ;  // internal exon ternimal part

    static Integer         k;            // order of the content MM
    static Integer         etorder;      // order of the exon terminating motif
    static Integer         etpseudocount;// pseudocount for the exon terminating motif
    static Integer         min_exon_length;
    static Integer         max_exon_length;
    static Integer         trans_init_window;

    static FramedPatMMGroup emiprobs;
    static FramedPatMMGroup highemiprobs;
    static FramedPatMMGroup mediumemiprobs;
    static FramedPatMMGroup lowemiprobs;
    
    static FramedPatMMGroup *GCemiprobs;
    
    static FramedPatMMGroup *highGCemiprobs;
    static FramedPatMMGroup *mediumGCemiprobs;
    static FramedPatMMGroup *lowGCemiprobs;
    
    static vector<Double>   initemiprobs[3];

    static vector<Double>   highinitemiprobs[3];
    static vector<Double>   mediuminitemiprobs[3];
    static vector<Double>   lowinitemiprobs[3];

    static vector<Double>   **GCinitemiprobs;

    static vector<Double>   **highGCinitemiprobs;
    static vector<Double>   **mediumGCinitemiprobs;
    static vector<Double>   **lowGCinitemiprobs;

    static vector<Double>   etemiprobs[3];

    static vector<Double>   highetemiprobs[3];
    static vector<Double>   mediumetemiprobs[3];
    static vector<Double>   lowetemiprobs[3];

    static vector<Double>   **GCetemiprobs;

    static vector<Double>   **highGCetemiprobs;
    static vector<Double>   **mediumGCetemiprobs;
    static vector<Double>   **lowGCetemiprobs;


    static vector<Integer>  numExonsOfType;
    static vector<Integer>  numHugeExonsOfType;  // number of exons exceeding the maximal length 
                                                 // modelled by the length distribution
    static vector<Integer>  highnumExonsOfType;
    static vector<Integer>  highnumHugeExonsOfType;  // number of exons exceeding the maximal length
                                                     // modelled by the length distribution

    static vector<Integer>  mediumnumExonsOfType;
    static vector<Integer>  mediumnumHugeExonsOfType;  // number of exons exceeding the maximal length
                                                     // modelled by the length distribution

    static vector<Integer>  lownumExonsOfType;
    static vector<Integer>  lownumHugeExonsOfType;  // number of exons exceeding the maximal length
                                                         // modelled by the length distribution

    static vector<Double> lenDistSingle;   // Length distribution of Single exons (length of biol. exon)
    static vector<Double> lenDistInitial;  // Length distribution of Initial exons (length of biol. exon)
    static vector<Double> lenDistInternal; // Length distribution of Internal exons (length of biol. exon)
    static vector<Double> lenDistTerminal; // Length distribution of Terminal exons (length of biol. exon)
    
    static vector<Double> highlenDistSingle;   // Length distribution of Single exons (length of biol. exon)
    static vector<Double> highlenDistInitial;  // Length distribution of Initial exons (length of biol. exon)
    static vector<Double> highlenDistInternal; // Length distribution of Internal exons (length of biol. exon)
    static vector<Double> highlenDistTerminal; // Length distribution of Terminal exons (length of biol. exon)
    
    static vector<Double> mediumlenDistSingle;   // Length distribution of Single exons (length of biol. exon)
    static vector<Double> mediumlenDistInitial;  // Length distribution of Initial exons (length of biol. exon)
    static vector<Double> mediumlenDistInternal; // Length distribution of Internal exons (length of biol. exon)
    static vector<Double> mediumlenDistTerminal; // Length distribution of Terminal exons (length of biol. exon)

    static vector<Double> lowlenDistSingle;   // Length distribution of Single exons (length of biol. exon)
    static vector<Double> lowlenDistInitial;  // Length distribution of Initial exons (length of biol. exon)
    static vector<Double> lowlenDistInternal; // Length distribution of Internal exons (length of biol. exon)
    static vector<Double> lowlenDistTerminal; // Length distribution of Terminal exons (length of biol. exon)
    
    static Integer numSingle, numInitial, numInternal, numTerminal; 
    
    static Integer highnumSingle, highnumInitial, highnumInternal, highnumTerminal;
    static Integer mediumnumSingle, mediumnumInitial, mediumnumInternal, mediumnumTerminal;
    static Integer lownumSingle, lownumInitial, lownumInternal, lownumTerminal;
    
    static Integer numHugeSingle, numHugeInitial, numHugeInternal, numHugeTerminal;

    static Integer highnumHugeSingle, highnumHugeInitial, highnumHugeInternal, highnumHugeTerminal;
    static Integer mediumnumHugeSingle, mediumnumHugeInitial, mediumnumHugeInternal, mediumnumHugeTerminal;
    static Integer lownumHugeSingle, lownumHugeInitial, lownumHugeInternal, lownumHugeTerminal;

    static Matrix<vector<Double> > Pls;
    
    static Matrix<vector<Double> > highPls;
    static Matrix<vector<Double> > mediumPls;
    static Matrix<vector<Double> > lowPls;
    
    static Matrix<vector<Double> >* GCPls;    // array with one matrix per GC content class
    
    static Matrix<vector<Double> >* highGCPls;
    static Matrix<vector<Double> >* mediumGCPls;
    static Matrix<vector<Double> >* lowGCPls;
    
    static Integer        exoncount;
    static Boolean        hasLenDist;
    static Integer        gesbasen[3];
    static Double         patpseudo;         // pseudocount for patterns in sequence
    static Integer        exonLenD;          // number of exons of length <= d
    static Integer        minPatSum;         // for the decision to shorten the emission pattern
    static vector<Integer> lenCountSingle;   // Length count of Single exons (length of biol. exon)

    static vector<Integer> highlenCountSingle;   // Length count of Single exons (length of biol. exon)
    static vector<Integer> mediumlenCountSingle;   // Length count of Single exons (length of biol. exon)
    static vector<Integer> lowlenCountSingle;   // Length count of Single exons (length of biol. exon)

    static vector<Integer> lenCountInitial;  // Length count of Initial exons (length of biol. exon)

    static vector<Integer> highlenCountInitial;  // Length count of Initial exons (length of biol. exon)
    static vector<Integer> mediumlenCountInitial;  // Length count of Initial exons (length of biol. exon)
    static vector<Integer> lowlenCountInitial;  // Length count of Initial exons (length of biol. exon)

    static vector<Integer> lenCountInternal; // Length count of Internal exons (length of biol. exon)

    static vector<Integer> highlenCountInternal; // Length count of Internal exons (length of biol. exon)
    static vector<Integer> mediumlenCountInternal; // Length count of Internal exons (length of biol. exon)
    static vector<Integer> lowlenCountInternal; // Length count of Internal exons (length of biol. exon)

    static vector<Integer> lenCountTerminal; // Length count of Terminal exons (length of biol. exon)

    static vector<Integer> highlenCountTerminal; // Length count of Terminal exons (length of biol. exon)
    static vector<Integer> mediumlenCountTerminal; // Length count of Terminal exons (length of biol. exon)
    static vector<Integer> lowlenCountTerminal; // Length count of Terminal exons (length of biol. exon)

    static double         slope_of_bandwidth;// for smoothing
    static Integer        minwindowcount;    // see class Smooth in commontrain.hh
    static Motif          *transInitMotif;   // weight matrix before the translation initiation
    static Motif          *GCtransInitMotif; // array for each GC content class
    static Integer        tis_motif_memory;  // order of the trans init motif
    static Integer        tis_motif_radius;  // radius for the smoothing of the trans init motif
    static Motif          **etMotif;         // weight matrices before the donor splice site (3 frames)
    static Motif          ***GCetMotif;      // array with motifs for each GC content class
    static int            numModels;
    static Double         *modelStartProbs;
    static int            ilend;
    static OpenReadingFrame *orf;
    static int            ochrecount, ambercount, opalcount; // frequencies of the 3 stop codons
    static bool           initAlgorithmsCalled, haveORF;
    static int            lastParIndex; // GC-index of current parameter set   
    static int            verbosity;
};


#endif  //  _EXONMODEL_HH
