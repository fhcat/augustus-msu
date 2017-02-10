 /*****************************************************************************\
 * Filename : utrmodel.hh
 * Author   : Mario Stanke
 * Description: Untranslated Region Model Header File
 *
 *
 * Date       |   Author        |  Changes
 *------------|-----------------|----------------------------------------
 * 21.09.2005 | Mario Stanke    | creation of the file
 * 27.03.2006 | Mario Stanke    | introduced UTR intron
 * 05.04.2007 | Mario Stanke    | distance distribution of tata box to tss
 \******************************************************************************/

#ifndef _UTRMODEL_HH
#define _UTRMODEL_HH

#include "statemodel.hh"


/**
 * The utr model class.
 *
 * author Mario Stanke
 */
 
class UtrModel : public StateModel {
public:
    UtrModel();
    ~UtrModel();

    StateType getStateType( ) const {
	return utype;
    }

  /**
   * Build all needed probabilities from the given
   * gene set.
   *
   * @param   annoseq A single linked list annotated sequences
   */
  void buildModel( const AnnoSequence* annoseq, int parIndex );
  void registerPars( Parameters* parameters);
  void processStates( const Gene* gene );
  void process5SingleExon( State* exon, bool withLen=true );
  void process5InitialExon(  State* exon, bool withLen=true );
  void process5InternalExon(  State* exon);
  void process5TerminalExon( State* exon);
  void process5Intron( int begin, int end);
  void process3SingleExon(  State* exon, bool withLen=true );
  void process3InitialExon(  State* exon, bool withLen=true );
  void process3InternalExon(  State* exon);
  void process3TerminalExon(  State* exon, bool withLen);
  void process3Intron( int begin, int end);
  
  
    /**
     * Print the wanted probabilities in files.
     * The names of the files should be given with
     * the Properties object /UtrModel/outfile
     */
    void printProbabilities   ( int zusNumber, BaseCount *bc, const char* suffix = NULL );
    void initAlgorithms       ( Matrix<Double>&, int);
    void viterbiForwardAndSampling(ViterbiMatrixType&, ViterbiMatrixType&, int, int, 
				   AlgorithmVariant, OptionListItem&) const;
    Double emiProbUnderModel  (int begin, int end) const;
    Double endPartEmiProb     (int begin, int end, int endOfBioExon) const;
    Double notEndPartEmiProb  (int begin, int end, int endOfBioExon, Feature *exonparts) const;
    void getEndPositions      ( int end, int &beginOfEndPart, int &endOfBioExon) const;
    Double tssupSeqProb       ( int left, int right, bool reverse) const;  
    Double tssProb            ( int left) const;
    void computeTtsProbs    ( );
    static void init();
    static void resetPars(){
	if (utrcount == 0)
	    return;
	if (!haveSnippetProbs)
	    initSnippetProbs();
	initAlgorithmsCalled = false;
    }
    static void readProbabilities(int zusNumber);
    static void readAllParameters();
    static void storeGCPars(int idx);
    static void resetModelCount(){utrcount = 0;};
    static void setTtsSpacing(int spacing){ ttsSpacing = spacing; };
    
private:
  Double seqProb            ( int left, int right, bool reverse, int type) const;
  void computeLengthDistributions( );
  static void fillTailsOfLengthDistributions( );
  void process5InitSequence( const char* start, const char* end);
  void process5Sequence( const char* start, const char* end);
  void process3Sequence( const char* start, const char* end);
  void processTssupSequence( const char* start, const char* end);
  /**
   *
   */
  void buildProbabilities ( const AnnoSequence* annoseq );
  void buildTSSModel( const AnnoSequence* annoseq );
  void buildTTSModel( const AnnoSequence* annoseq );
  int findTATA(const char* seq, int maxpos, bool reverseComplement=false) const;
  void processTSS(const char* start);
  //  void buildLenDist       ( const AnnoSequence* annoseq );
  //  void storeUtrLengths ( const AnnoSequence* annoseq);
  void initCountVars      ( );
  Double longIntronProb(int internalBegin, int internalEnd) const;
  static void initSnippetProbs();

private:
  StateType              utype;
  Integer                gweight;
  static Integer         utrcount;
  static Integer         utrcountHigh;
  static Integer         utrcountMedium;
  static Integer         utrcountLow;
  
  static vector<Integer> utr5_emicount;
  static vector<Integer> utr5_emicountHigh;
  static vector<Integer> utr5_emicountMedium;
  static vector<Integer> utr5_emicountLow;
  
  static vector<Integer> utr5init_emicount;
  static vector<Integer> utr5init_emicountHigh;
  static vector<Integer> utr5init_emicountMedium;
  static vector<Integer> utr5init_emicountLow;
  
  static vector<Integer> utr3_emicount;
  static vector<Integer> utr3_emicountHigh;
  static vector<Integer> utr3_emicountMedium;
  static vector<Integer> utr3_emicountLow;
  
  static Double          utr_patpseudo;
  static Double          utr_patpseudoHigh;
  static Double          utr_patpseudoMedium;
  static Double          utr_patpseudoLow;
  
  static PatMMGroup      utr5_emiprobs;
  static PatMMGroup      utr5_emiprobsHigh;
  static PatMMGroup      utr5_emiprobsMedium;
  static PatMMGroup      utr5_emiprobsLow;
  
  static PatMMGroup      *GCutr5_emiprobs;
  static PatMMGroup      *GCutr5_emiprobsHigh;
  static PatMMGroup      *GCutr5_emiprobsMedium;
  static PatMMGroup      *GCutr5_emiprobsLow;
  
  static PatMMGroup      utr5init_emiprobs;
  static PatMMGroup      utr5init_emiprobsHigh;
  static PatMMGroup      utr5init_emiprobsMedium;
  static PatMMGroup      utr5init_emiprobsLow;
  
  static PatMMGroup      *GCutr5init_emiprobs;
  static PatMMGroup      *GCutr5init_emiprobsHigh;
  static PatMMGroup      *GCutr5init_emiprobsMedium;
  static PatMMGroup      *GCutr5init_emiprobsLow;
  
  static PatMMGroup      utr3_emiprobs;
  static PatMMGroup      utr3_emiprobsHigh;
  static PatMMGroup      utr3_emiprobsMedium;
  static PatMMGroup      utr3_emiprobsLow;
  
  static PatMMGroup      *GCutr3_emiprobs;
  static PatMMGroup      *GCutr3_emiprobsHigh;
  static PatMMGroup      *GCutr3_emiprobsMedium;
  static PatMMGroup      *GCutr3_emiprobsLow;
  
  static Integer         utr5init_gesbasen;
  static Integer         utr5init_gesbasenHigh;
  static Integer         utr5init_gesbasenMedium;
  static Integer         utr5init_gesbasenLow;
  
  static Integer         utr5_gesbasen;
  static Integer         utr5_gesbasenHigh;
  static Integer         utr5_gesbasenMedium;
  static Integer         utr5_gesbasenLow;
  
  static Integer         utr3_gesbasen;
  static Integer         utr3_gesbasenHigh;
  static Integer         utr3_gesbasenMedium;
  static Integer         utr3_gesbasenLow;
  
  static Integer         k;
  static double          utr5patternweight; // old way: this is applied AFTER reading from the parameters file
  static double          utr5patternweightHigh; // old way: this is applied AFTER reading from the parameters file
  static double          utr5patternweightMedium; // old way: this is applied AFTER reading from the parameters file
  static double          utr5patternweightLow; // old way: this is applied AFTER reading from the parameters file
   
  static double          utr3patternweight; // old way: this is applied AFTER reading from the parameters file
  static double          utr3patternweightHigh; // old way: this is applied AFTER reading from the parameters file
  static double          utr3patternweightMedium; // old way: this is applied AFTER reading from the parameters file
  static double          utr3patternweightLow; // old way: this is applied AFTER reading from the parameters file
  
  static double          utr5prepatternweight; // for computing a mixture directly after HMM training and BEFORE writing down to parameter file
  static double          utr5prepatternweightHigh; // for computing a mixture directly after HMM training and BEFORE writing down to parameter file
  static double          utr5prepatternweightMedium; // for computing a mixture directly after HMM training and BEFORE writing down to parameter file
  static double          utr5prepatternweightLow; // for computing a mixture directly after HMM training and BEFORE writing down to parameter file
  
  static double          utr3prepatternweight; // for computing a mixture directly after HMM training and BEFORE writing down to parameter file
  static double          utr3prepatternweightHigh; // for computing a mixture directly after HMM training and BEFORE writing down to parameter file
  static double          utr3prepatternweightMedium; // for computing a mixture directly after HMM training and BEFORE writing down to parameter file
  static double          utr3prepatternweightLow; // for computing a mixture directly after HMM training and BEFORE writing down to parameter file
  
  static vector<Integer> tssup_emicount;
  static Double          tssup_patpseudo;
  static vector<Double>  tssup_emiprobs;
  static vector<Double>  *GCtssup_emiprobs;
  static Integer         tssup_gesbasen;
  static Integer         tssup_k;
  
  static vector<Integer> lenCount5Single;       // Length count of single exons
  static vector<Integer> lenCount5SingleHigh;       // Length count of single exons
  static vector<Integer> lenCount5SingleMedium;       // Length count of single exons
  static vector<Integer> lenCount5SingleLow;       // Length count of single exons
  
  static vector<Integer> lenCount5Initial;      // Length count of initial exons
  static vector<Integer> lenCount5InitialHigh;      // Length count of initial exons
  static vector<Integer> lenCount5InitialMedium;      // Length count of initial exons
  static vector<Integer> lenCount5InitialLow;      // Length count of initial exons
  
  static vector<Integer> lenCount5Internal;     // Length count of internal exons
  static vector<Integer> lenCount5InternalHigh;     // Length count of internal exons
  static vector<Integer> lenCount5InternalMedium;     // Length count of internal exons
  static vector<Integer> lenCount5InternalLow;     // Length count of internal exons
  
  static vector<Integer> lenCount5Terminal;     // Length count of terminal exons
  static vector<Integer> lenCount5TerminalHigh;     // Length count of terminal exons
  static vector<Integer> lenCount5TerminalMedium;     // Length count of terminal exons
  static vector<Integer> lenCount5TerminalLow;     // Length count of terminal exons
  
  static vector<Double>  lenDist5Single;        // Length distribution of single exons
  static vector<Double>  lenDist5SingleHigh;        // Length distribution of single exons
  static vector<Double>  lenDist5SingleMedium;        // Length distribution of single exons
  static vector<Double>  lenDist5SingleLow;        // Length distribution of single exons
  
  
  static vector<Double>  lenDist5Initial;       // Length distribution of initial exons
  static vector<Double>  lenDist5InitialHigh;       // Length distribution of initial exons
  static vector<Double>  lenDist5InitialMedium;       // Length distribution of initial exons
  static vector<Double>  lenDist5InitialLow;       // Length distribution of initial exons
  
  static vector<Double>  lenDist5Internal;      // Length distribution of internal exons
  static vector<Double>  lenDist5InternalHigh;      // Length distribution of internal exons
  static vector<Double>  lenDist5InternalMedium;      // Length distribution of internal exons
  static vector<Double>  lenDist5InternalLow;      // Length distribution of internal exons
  
  
  static vector<Double>  lenDist5Terminal;      // Length distribution of terminal exons
  static vector<Double>  lenDist5TerminalHigh;      // Length distribution of terminal exons
  static vector<Double>  lenDist5TerminalMedium;      // Length distribution of terminal exons
  static vector<Double>  lenDist5TerminalLow;      // Length distribution of terminal exons
  
  static vector<Double>  tailLenDist5Single;    // Tail probabilities of the length distribution of single exons
  static vector<Double>  tailLenDist5SingleHigh;    // Tail probabilities of the length distribution of single exons
  static vector<Double>  tailLenDist5SingleMedium;    // Tail probabilities of the length distribution of single exons
  static vector<Double>  tailLenDist5SingleLow;    // Tail probabilities of the length distribution of single exons
  
  static vector<Integer> lenCount3Single;       // Length count of single exons
  static vector<Integer> lenCount3SingleHigh;       // Length count of single exons
  static vector<Integer> lenCount3SingleMedium;       // Length count of single exons
  static vector<Integer> lenCount3SingleLow;       // Length count of single exons
  
  static vector<Integer> lenCount3Initial;      // Length count of initial exons
  static vector<Integer> lenCount3InitialHigh;      // Length count of initial exons
  static vector<Integer> lenCount3InitialMedium;      // Length count of initial exons
  static vector<Integer> lenCount3InitialLow;      // Length count of initial exons
  
  static vector<Integer> lenCount3Internal;     // Length count of internal exons
  
  static vector<Integer> lenCount3InternalHigh;     // Length count of internal exons
  static vector<Integer> lenCount3InternalMedium;     // Length count of internal exons
  static vector<Integer> lenCount3InternalLow;     // Length count of internal exons
  
  static vector<Integer> lenCount3Terminal;     // Length count of terminal exons
  static vector<Integer> lenCount3TerminalHigh;     // Length count of terminal exons
  static vector<Integer> lenCount3TerminalMedium;     // Length count of terminal exons
  static vector<Integer> lenCount3TerminalLow;     // Length count of terminal exons
  
  static vector<Double>  lenDist3Single;        // Length distribution of single exons
  static vector<Double>  lenDist3SingleHigh;        // Length distribution of single exons
  static vector<Double>  lenDist3SingleMedium;        // Length distribution of single exons
  static vector<Double>  lenDist3SingleLow;        // Length distribution of single exons
  
  static vector<Double>  lenDist3Initial;       // Length distribution of initial exons
  static vector<Double>  lenDist3InitialHigh;       // Length distribution of initial exons
  static vector<Double>  lenDist3InitialMedium;       // Length distribution of initial exons
  static vector<Double>  lenDist3InitialLow;       // Length distribution of initial exons
  
  static vector<Double>  lenDist3Internal;      // Length distribution of internal exons
  static vector<Double>  lenDist3InternalHigh;      // Length distribution of internal exons
  static vector<Double>  lenDist3InternalMedium;      // Length distribution of internal exons
  static vector<Double>  lenDist3InternalLow;      // Length distribution of internal exons
  
  static vector<Double>  lenDist3Terminal;      // Length distribution of terminal exons
  static vector<Double>  lenDist3TerminalHigh;      // Length distribution of terminal exons
  static vector<Double>  lenDist3TerminalMedium;      // Length distribution of terminal exons
  static vector<Double>  lenDist3TerminalLow;      // Length distribution of terminal exons
  
  
  static vector<Double>  tailLenDist3Single;    // Tail probabilities of the length distribution of single exons
  static vector<Double>  tailLenDist3SingleHigh;    // Tail probabilities of the length distribution of single exons
  static vector<Double>  tailLenDist3SingleMedium;    // Tail probabilities of the length distribution of single exons
  static vector<Double>  tailLenDist3SingleLow;    // Tail probabilities of the length distribution of single exons
  
  
  static vector<Double>  tssProbsPlus;          // to store tss probabilities
  static vector<Double>  tssProbsMinus;         // to store tss probabilities
  
  static Integer         num5Single, num5SingleHigh, num5SingleMedium, num5SingleLow, num5Initial, num5InitialHigh, num5InitialMedium, num5InitialLow, num5Internal, num5InternalHigh, num5InternalMedium, num5InternalLow, num5Terminal, num5TerminalHigh, num5TerminalMedium, num5TerminalLow, num5Introns;
  static Integer         numHuge5Single, numHuge5SingleHigh, numHuge5SingleMedium, numHuge5SingleLow, numHuge5Initial, numHuge5InitialHigh, numHuge5InitialMedium, numHuge5InitialLow, numHuge5Internal, numHuge5InternalHigh, numHuge5InternalMedium, numHuge5InternalLow, numHuge5Terminal, numHuge5TerminalHigh, numHuge5TerminalMedium, numHuge5TerminalLow; 
  static Integer         num3Single, num3SingleHigh, num3SingleMedium, num3SingleLow, num3Initial, num3InitialHigh, num3InitialMedium, num3InitialLow, num3Internal, num3InternalHigh, num3InternalMedium, num3InternalLow, num3Terminal, num3TerminalHigh, num3TerminalMedium, num3TerminalLow, num3Introns;
  static Integer         numHuge3Single, numHuge3SingleHigh, numHuge3SingleMedium, numHuge3SingleLow, numHuge3Initial, numHuge3InitialHigh, numHuge3InitialMedium, numHuge3InitialLow, numHuge3Internal, numHuge3InternalHigh, numHuge3InternalMedium, numHuge3InternalLow, numHuge3Terminal, numHuge3TerminalHigh, numHuge3TerminalMedium, numHuge3TerminalLow; 
  static Integer         exonLenD, exonLenDHigh, exonLenDMedium, exonLenDLow;            // use detailed length distribution up to this number
  static Integer         max_exon_length, max_exon_lengthHigh, max_exon_lengthMedium, max_exon_lengthLow;
  static Integer         max3singlelength, max3singlelengthHigh, max3singlelengthMedium, max3singlelengthLow;
  static Integer         max3termlength, max3termlengthHigh, max3termlengthMedium, max3termlengthLow;
  static double          slope_of_bandwidth;  // for smoothing
  static Integer         minwindowcount;      // see class Smooth in commontrain.hh
  static Boolean         hasLenDist;
  static Integer         tss_start;
  static Integer         tss_end;
  static Integer         tata_start;
  static Integer         tata_end;
  static Integer         tata_pseudocount;
  static Integer         d_tss_tata_min;
  static Integer         d_tss_tata_max;
  static Motif           *tssMotif;           // motif of the transcription start site of tata-less promotors
  static Motif           *GCtssMotif;
  static Motif           *ttsMotif;           // motif of the transcription termination site (downstream of polyA signal)
  static Motif           *GCttsMotif;
  static Motif           *tssMotifTATA;       // motif of the transcription start site of tata promotors
  static Motif           *GCtssMotifTATA;
  static Motif           *tataMotif;          // motif of the tata box (if existent)
  static Motif           *GCtataMotif;
  // UTR intron related member variables
  static vector<Integer> intron_emicount;
  //static vector<Double>  intron_emiprobs;
  //static Integer         intron_k;            // order of the markov chain
  static SnippetProbs    *initSnippetProbs5, *snippetProbs5, *rInitSnippetProbs5, *rSnippetProbs5, *rSnippetProbs3, *intronSnippetProbs;

  static SnippetProbs    *initSnippetProbs5High;
  static SnippetProbs    *initSnippetProbs5Medium;
  static SnippetProbs    *initSnippetProbs5Low;

  static SnippetProbs    *snippetProbs5High;
  static SnippetProbs    *snippetProbs5Medium;
  static SnippetProbs    *snippetProbs5Low;

  static SnippetProbs    *rSnippetProbs5High;
  static SnippetProbs    *rSnippetProbs5Medium;
  static SnippetProbs    *rSnippetProbs5Low;

  static SnippetProbs    *rInitSnippetProbs5High;
  static SnippetProbs    *rInitSnippetProbs5Medium;
  static SnippetProbs    *rInitSnippetProbs5Low;

  static SnippetProbs    *rSnippetProbs3High;
  static SnippetProbs    *rSnippetProbs3Medium;
  static SnippetProbs    *rSnippetProbs3Low;

  static bool            initAlgorithmsCalled, haveSnippetProbs;
  static vector<Integer> aataaa_count;
  static vector<Double>  aataaa_probs;
  static int             aataaa_boxlen;
  static string          polyasig_consensus;
  static int             d_polya_cleavage_min;
  static int             d_polya_cleavage_max;
  static double          prob_polya;
  static int             tts_motif_memory;
  static double pUtr5Intron, pUtr3Intron, prUtr5Intron, prUtr3Intron;
  static Double          *ttsProbPlus, *ttsProbMinus;
  static vector<Integer> distCountTata;
  static int             lastParIndex;
  static int             verbosity;
  static int             ttsSpacing; // without hints allow 3' end only every ttsSpacing bases for speed
};

class UtrModelError : public ProjectError {
public:
    UtrModelError(string msg) : ProjectError(msg) {}
};

#endif    //  _UTRMODEL_HH
