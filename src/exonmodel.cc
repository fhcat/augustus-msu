/**********************************************************************
 * file:    exonmodel.cc
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  coding exons only, UTR exons are modelled in utrmodel.cc
 * authors: Mario Stanke, mario@gobics.de
 *
 * date    |   author      |  changes
 * --------|---------------|------------------------------------------
 * 06.05.02| Mario Stanke  | debugging and testing
 * 18.05.02| Mario Stanke  | added length distributions
 * 31.05.02| Mario Stanke  | implementation of viterbi algorithm
 * 06.10.02| Mario Stanke  | making length distr independent of content
 * 18.11.02| Mario Stanke  | internal exon terminal emission part
 * 30.12.02| Mario Stanke  | moving splice site models to intron
 * 25.03.03| Mario Stanke  | shadow exons on reverse strand
 * 02.07.04| Mario Stanke  | malus for all exons, not only for those without a bonus
 * 19.08.04| Mario Stanke  | corrected rare bug with exon state at end of dna
 * 23.12.04| Mario Stanke  | bug fixed: wrong reading frame in call to getExonListInRange
 * 08.02.05| Mario Stanke  | bug fixed: exon hints were missed for reverse terminal or reverse single exons.
 * 22.08.05| Mario Stanke  | training on genes with incomplete 3' end, incomplete exon not used for lendist
 * 09.09.08| Mario Stanke  | introduce fuzzyness in dss and ass hints
 **********************************************************************/

#include "exonmodel.hh"

// project includes
#include "motif.hh"
#include "pp_scoring.hh"
#include "properties.hh"
#include "projectio.hh"
#include "extrinsicinfo.hh"

// standard C/C++ includes
#include <fstream>

/*
 * Initialisation of static data members
 */

vector<Integer> ExonModel::patterncount[3];      // {0,1,2}x{acgt}^(k+1), the reading frame is 
                                                 // the positiion of the emitted (last) nucleotide
vector<Integer> ExonModel::patternhighcount[3];
vector<Integer> ExonModel::patternmediumcount[3];
vector<Integer> ExonModel::patternlowcount[3];
                                                 
vector<Integer> ExonModel::initpatterncount[3];

FramedPatMMGroup ExonModel::highemiprobs("exon emiprob old");
FramedPatMMGroup ExonModel::mediumemiprobs("exon emiprob old");
FramedPatMMGroup ExonModel::lowemiprobs("exon emiprob old");

vector<Integer> ExonModel::highinitpatterncount[3];
vector<Integer> ExonModel::mediuminitpatterncount[3];
vector<Integer> ExonModel::lowinitpatterncount[3];

vector<Integer> ExonModel::etpatterncount[3];

vector<Integer> ExonModel::highetpatterncount[3];
vector<Integer> ExonModel::mediumetpatterncount[3];
vector<Integer> ExonModel::lowetpatterncount[3];

Integer ExonModel::k = 4;
Integer ExonModel::etorder = 1;
Integer ExonModel::etpseudocount = 3;
Integer ExonModel::min_exon_length = 1;
Integer ExonModel::max_exon_length = 12000;
Integer ExonModel::trans_init_window = 12;
FramedPatMMGroup ExonModel::emiprobs("exon emiprob old");


FramedPatMMGroup* ExonModel::GCemiprobs = NULL;

FramedPatMMGroup* ExonModel::highGCemiprobs = NULL;
FramedPatMMGroup* ExonModel::mediumGCemiprobs = NULL;
FramedPatMMGroup* ExonModel::lowGCemiprobs = NULL;

vector<Double> ExonModel::initemiprobs[3];

vector<Double> ExonModel::highinitemiprobs[3];
vector<Double> ExonModel::mediuminitemiprobs[3];
vector<Double> ExonModel::lowinitemiprobs[3];

vector<Double> **ExonModel::GCinitemiprobs = NULL;

vector<Double> **ExonModel::highGCinitemiprobs = NULL;
vector<Double> **ExonModel::mediumGCinitemiprobs = NULL;
vector<Double> **ExonModel::lowGCinitemiprobs = NULL;

vector<Double> ExonModel::etemiprobs[3];

vector<Double> ExonModel::highetemiprobs[3];
vector<Double> ExonModel::mediumetemiprobs[3];
vector<Double> ExonModel::lowetemiprobs[3];


vector<Double> **ExonModel::GCetemiprobs = NULL;

vector<Double> **ExonModel::highGCetemiprobs = NULL;
vector<Double> **ExonModel::mediumGCetemiprobs = NULL;
vector<Double> **ExonModel::lowGCetemiprobs = NULL;


vector<Integer> ExonModel::numExonsOfType;
vector<Integer> ExonModel::numHugeExonsOfType;  // number of exons exceeding the maximal length 
                                                    // modelled by the length distribution

vector<Integer> ExonModel::highnumExonsOfType;
vector<Integer> ExonModel::highnumHugeExonsOfType;  // number of exons exceeding the maximal length
                                                    // modelled by the length distribution

vector<Integer> ExonModel::mediumnumExonsOfType;
vector<Integer> ExonModel::mediumnumHugeExonsOfType;  // number of exons exceeding the maximal length
                                                    // modelled by the length distribution

vector<Integer> ExonModel::lownumExonsOfType;
vector<Integer> ExonModel::lownumHugeExonsOfType;  // number of exons exceeding the maximal length
                                                    // modelled by the length distribution

vector<Double> ExonModel::lenDistSingle;   // Length distribution of Single exons 
                                               // (length of biol. ExonModel::exon)

vector<Double> ExonModel::highlenDistSingle;   // Length distribution of Single exons
                                               // (length of biol. ExonModel::exon)
vector<Double> ExonModel::mediumlenDistSingle;   // Length distribution of Single exons
                                               // (length of biol. ExonModel::exon)
vector<Double> ExonModel::lowlenDistSingle;   // Length distribution of Single exons
                                               // (length of biol. ExonModel::exon)

vector<Double> ExonModel::lenDistInitial;  // Length distribution of Initial exons 
                                               // (length of biol. exon)

vector<Double> ExonModel::highlenDistInitial;  // Length distribution of Initial exons
                                               // (length of biol. exon)
vector<Double> ExonModel::mediumlenDistInitial;  // Length distribution of Initial exons
                                               // (length of biol. exon)
vector<Double> ExonModel::lowlenDistInitial;  // Length distribution of Initial exons
                                               // (length of biol. exon)

vector<Double> ExonModel::lenDistInternal; // Length distribution of Internal exons 
                                               // (length of biol. exon)

vector<Double> ExonModel::highlenDistInternal; // Length distribution of Internal exons
                                               // (length of biol. exon)
vector<Double> ExonModel::mediumlenDistInternal; // Length distribution of Internal exons
                                               // (length of biol. exon)
vector<Double> ExonModel::lowlenDistInternal; // Length distribution of Internal exons
                                               // (length of biol. exon)

vector<Double> ExonModel::lenDistTerminal; // Length distribution of Terminal exons 
                                               // (length of biol. exon)

vector<Double> ExonModel::highlenDistTerminal; // Length distribution of Terminal exons
                                               // (length of biol. exon)
vector<Double> ExonModel::mediumlenDistTerminal; // Length distribution of Terminal exons
                                               // (length of biol. exon)
vector<Double> ExonModel::lowlenDistTerminal; // Length distribution of Terminal exons
                                               // (length of biol. exon)

Motif*   ExonModel::transInitMotif = NULL;      // weight matrix before the translation initiation
Motif*   ExonModel::GCtransInitMotif = NULL;    // weight matrix before the translation initiation
Motif**  ExonModel::etMotif = NULL;             // weight matrices before the donor splice site (3 frames)
Motif*** ExonModel::GCetMotif = NULL;           // array of above
Integer  ExonModel::numSingle=0, ExonModel::numInitial=0, ExonModel::numInternal=0, 
         ExonModel::numTerminal=0;
Integer  ExonModel::highnumSingle=0, ExonModel::highnumInitial=0, ExonModel::highnumInternal=0,
         ExonModel::highnumTerminal=0;
Integer  ExonModel::mediumnumSingle=0, ExonModel::mediumnumInitial=0, ExonModel::mediumnumInternal=0,
         ExonModel::mediumnumTerminal=0;
Integer  ExonModel::lownumSingle=0, ExonModel::lownumInitial=0, ExonModel::lownumInternal=0,
         ExonModel::lownumTerminal=0;
Integer  ExonModel::numHugeSingle=0, ExonModel::numHugeInitial=0, 
         ExonModel::numHugeInternal=0, ExonModel::numHugeTerminal=0;

Integer  ExonModel::highnumHugeSingle=0, ExonModel::highnumHugeInitial=0,
         ExonModel::highnumHugeInternal=0, ExonModel::highnumHugeTerminal=0;
Integer  ExonModel::mediumnumHugeSingle=0, ExonModel::mediumnumHugeInitial=0,
         ExonModel::mediumnumHugeInternal=0, ExonModel::mediumnumHugeTerminal=0;
Integer  ExonModel::lownumHugeSingle=0, ExonModel::lownumHugeInitial=0,
         ExonModel::lownumHugeInternal=0, ExonModel::lownumHugeTerminal=0;

Matrix<vector<Double> > ExonModel::Pls;

Matrix<vector<Double> > ExonModel::highPls;
Matrix<vector<Double> > ExonModel::mediumPls;
Matrix<vector<Double> > ExonModel::lowPls;

Matrix<vector<Double> >* ExonModel::GCPls = NULL;

Matrix<vector<Double> >* ExonModel::highGCPls = NULL;
Matrix<vector<Double> >* ExonModel::mediumGCPls = NULL;
Matrix<vector<Double> >* ExonModel::lowGCPls = NULL;

Integer ExonModel::exoncount = 0;
Boolean ExonModel::hasLenDist = false;
//Boolean ExonModel::hasAAdep = false;      // not in use right now
Integer ExonModel::gesbasen[3]  = { 0, 0, 0 };
Double  ExonModel::patpseudo = 1;         // pseudocount for patterns in sequence
Integer ExonModel::exonLenD = 4000;       // use detailled length distribution up to this number
Integer ExonModel::minPatSum = 0;         // for the decision to shorten the emission pattern
vector<Integer> ExonModel::lenCountSingle;   // length count of Single exons (length of biol. exon)

vector<Integer> ExonModel::highlenCountSingle;   // length count of Single exons (length of biol. exon)
vector<Integer> ExonModel::mediumlenCountSingle;   // length count of Single exons (length of biol. exon)
vector<Integer> ExonModel::lowlenCountSingle;   // length count of Single exons (length of biol. exon)

vector<Integer> ExonModel::lenCountInitial;  // length count of Initial exons (length of biol. exon)

vector<Integer> ExonModel::highlenCountInitial;  // length count of Initial exons (length of biol. exon)
vector<Integer> ExonModel::mediumlenCountInitial;  // length count of Initial exons (length of biol. exon)
vector<Integer> ExonModel::lowlenCountInitial;  // length count of Initial exons (length of biol. exon)

vector<Integer> ExonModel::lenCountInternal; // length count of Internal exons (length of biol. exon)

vector<Integer> ExonModel::highlenCountInternal; // length count of Internal exons (length of biol. exon)
vector<Integer> ExonModel::mediumlenCountInternal; // length count of Internal exons (length of biol. exon)
vector<Integer> ExonModel::lowlenCountInternal; // length count of Internal exons (length of biol. exon)

vector<Integer> ExonModel::lenCountTerminal; // length count of Terminal exons (length of biol. exon)

vector<Integer> ExonModel::highlenCountTerminal; // length count of Terminal exons (length of biol. exon)
vector<Integer> ExonModel::mediumlenCountTerminal; // length count of Terminal exons (length of biol. exon)
vector<Integer> ExonModel::lowlenCountTerminal; // length count of Terminal exons (length of biol. exon)

double ExonModel::slope_of_bandwidth = 0.1;     // for smoothing
Integer ExonModel::minwindowcount = 5;           // see class Smooth in commontrain.hh
Integer ExonModel::tis_motif_memory = 3;           // see class Smooth in commontrain.hh
Integer ExonModel::tis_motif_radius = 2;           // see class Smooth in commontrain.hh

int             ExonModel::numModels = 1;
//BaumWelch*      ExonModel::baumwelch = NULL;
Double*         ExonModel::modelStartProbs = NULL;
//AADependency    ExonModel::aadep = 0;
int             ExonModel::ilend = 550;
OpenReadingFrame* ExonModel::orf = NULL;
int             ExonModel::ochrecount = 0; // frequencies of the stop codons
int             ExonModel::ambercount = 0;
int             ExonModel::opalcount  = 0;
bool            ExonModel::initAlgorithmsCalled = false;
bool            ExonModel::haveORF = false;
int             ExonModel::lastParIndex = -1; // GC-index of current parameter set
int             ExonModel::verbosity;


/* --- OpenReadingFrame methods ------------------------------------ */

/*
 * constructor
 */
OpenReadingFrame::OpenReadingFrame(const char *dna, int _max_exon_length, int _n) :
    n(_n), max_exon_length(_max_exon_length)
{
    nearestStopForward.resize(n);
    nearestStopReverse.resize(n);
    int stopcodpos, i;
    stopcodpos=-1;
    for (i=0; i<=n - STOPCODON_LEN; i+=3) {
	if (isStopcodon(dna+i)) 
	    stopcodpos = i;
	nearestStopForward[i] = stopcodpos;
    }
    stopcodpos=-1;
    for (i=1; i<=n - STOPCODON_LEN; i+=3) {
	if (isStopcodon(dna+i))
	    stopcodpos = i;
	nearestStopForward[i] = stopcodpos;
    }
    stopcodpos=-1;
    for (i=2; i<=n - STOPCODON_LEN; i+=3) {
	if (isStopcodon(dna+i))
	    stopcodpos = i;
	nearestStopForward[i] = stopcodpos;
    }
    stopcodpos=-1;
    for (i=0; i<=n - STOPCODON_LEN; i+=3) {
	if (isRCStopcodon(dna+i))
	    stopcodpos = i;
	nearestStopReverse[i] = stopcodpos;
    }
    stopcodpos=-1;
    for (i=1; i<=n - STOPCODON_LEN; i+=3) {
	if (isRCStopcodon(dna+i))
	    stopcodpos = i;
	nearestStopReverse[i] = stopcodpos;
    }
    stopcodpos=-1;
    for (i=2; i<=n - STOPCODON_LEN; i+=3) {
	if (isRCStopcodon(dna+i))
	    stopcodpos = i;
	nearestStopReverse[i] = stopcodpos;
    }
    if (n>5) {
	// nearestStopForward[n - STOPCODON_LEN] = nearestStopForward[n - STOPCODON_LEN - 3];
	nearestStopForward[n - STOPCODON_LEN +1] = nearestStopForward[n - STOPCODON_LEN - 2];
	nearestStopForward[n - STOPCODON_LEN +2] = nearestStopForward[n - STOPCODON_LEN - 1];
	// nearestStopReverse[n - STOPCODON_LEN] = nearestStopReverse[n - STOPCODON_LEN - 3];
	nearestStopReverse[n - STOPCODON_LEN +1] = nearestStopReverse[n - STOPCODON_LEN - 2];
	nearestStopReverse[n - STOPCODON_LEN +2] = nearestStopReverse[n - STOPCODON_LEN - 1];
    } else {// TODO (unwichtig)
    }
}


/*
 * Get the position furthest to the left of base such that with respect to 
 * frame at position base there is no stop codon in the window from 
 * leftmostExonBegin to base. If forward is false the reverse complement is taken.
 * Initialize the tables of stopcodons.
 */
int OpenReadingFrame::leftmostExonBegin(int frame, int base, bool forward){
    int pos;
    if (forward) {
	if (frame==0 || frame==1)
	    pos=base-frame-3;
	else 
	    pos=base-frame;
    } else {
	if (frame==1 || frame==2)
	    pos=base+frame-5;
	else 
	    pos=base-2;
    }
    
    Integer leftmostbegin;
    if (pos >= n)
	pos -= 3*((pos-n+3)/3);
    // alternativ: pos = n-3 + (pos-n)%3

    if (pos >= 0) {
	leftmostbegin = forward? nearestStopForward[pos]+1 : nearestStopReverse[pos]+1;
    }
    else 
	leftmostbegin = 0;
   

    // this has the effect that no exon can be longer than max_exon_length
    // which is actually not a necessary condition but it is easier this way
    // The 10 is a security distance, because the biological exon can be longer than the inner sequence part 
    int max_allowed_len = max_exon_length - Constant::ass_upwindow_size - Constant::ass_start - ASS_MIDDLE - DSS_MIDDLE - Constant::dss_start;
    if (leftmostbegin < base-max_allowed_len)
      leftmostbegin = base-max_allowed_len;
    return leftmostbegin;
}

/*
 * containsInFrameStopcodon: not needed, can be implemented by GeneticCode
 */

bool OpenReadingFrame::isStopcodon( const char* dna ){
    if (*dna != tolower(*dna))
	throw ProjectError("Internal Error: upper case dna"); // TEMP: just for now
    // Oliver, 2009/06/24: commented this out to give translation table priority
    // over {ochre,amber,opal}prob values; this basically causes
    // stopcodons with prob=0 to appear neither prematurely nor as stopcodons
    // if( *dna == 't' ){
    // 	if( *(dna+1) == 'a' )
    // 	    return ((*(dna+2) == 'a' && ochre) || (*(dna+2) == 'g' && amber));
    // 	if( *(dna+1) == 'g' )
    // 	    return (*(dna+2) == 'a' && opal);
    // }
    // return false;

    // give the chosen translation table priority over {ochre,amber,opal}prob
    // 
    return GeneticCode::isStopcodon(dna);
}

bool OpenReadingFrame::isRCStopcodon( const char* dna ){
    return GeneticCode::isRCStopcodon(dna);
}


/* --- ExonModel methods ------------------------------------------- */

/*
 * constructor
 */
ExonModel::ExonModel() : gweight(1) {
    etype = toStateType( Properties::getProperty("/ExonModel/type", exoncount++) );
    win = stateReadingFrames[etype]; // cout << "etype:    " << etype << endl;
    switch (etype) {
	case singleGHigh: case singleGMedium: case singleGLow:
	case initial0High: case initial0Medium: case initial0Low:
	case initial1High: case initial1Medium: case initial1Low:
	case initial2High: case initial2Medium: case initial2Low:
	    // start codon at left end of exon
	    beginPartLen = STARTCODON_LEN + trans_init_window;
	    innerPartOffset = STARTCODON_LEN;
	    break;
	case rsingleGHigh: case rsingleGMedium: case rsingleGLow:
	case rterminal0High: case rterminal0Medium: case rterminal0Low: 
  case rterminal1High: case rterminal1Medium: case rterminal1Low:
  case rterminal2High: case rterminal2Medium: case rterminal2Low:  
	    // reverse stop codon at left end of exon
	    beginPartLen = innerPartOffset = STOPCODON_LEN;
	    break;
	default:   // splice site at left end of exon
	    beginPartLen = 0;
	    innerPartOffset = isOnFStrand(etype) ? Constant::ass_end : Constant::dss_start;
    }
    switch (etype) {
	case singleGHigh: case singleGMedium: case singleGLow: case terminalHigh: case terminalMedium: case terminalLow:
	    // stop codon at right end of exon
	    innerPartEndOffset = STOPCODON_LEN;
	    baseOffset = 0;
	    // endPartLen = STOPCODON_LEN;
	    break;
	case rsingleGHigh: case rsingleGMedium: case rsingleGLow: case rinitialHigh: case rinitialMedium: case rinitialLow:
	    // reverse start codon at right end of exon
	    innerPartEndOffset = STARTCODON_LEN;
	    baseOffset = -trans_init_window;
	    // endPartLen = STARTCODON_LEN + trans_init_window;
	    break;
	default:  // splice site at right end of exon
	    innerPartEndOffset = baseOffset = isOnFStrand(etype) ? Constant::dss_start : Constant::ass_end;
	    // endPartLen = 0;
    }
}

/*
 * destructor
 */
ExonModel::~ExonModel( ){
    // TODO: delete what is neccessary
    if( --exoncount == 0 )
	lastParIndex = -1;
}

/*
 * ===[ ExonModel initialisation of class variables ]======================
 */
void ExonModel::init() {
    /*
     * Default values for the parameters.
     */
    k = 4;
    patpseudo  = Double(1.0);   
    slope_of_bandwidth = 0.1;   // for smoothing
    minwindowcount=3;           // see class Smooth in commontrain.hh	
    try{
	verbosity = Properties::getIntProperty("/ExonModel/verbosity");
    }catch( ProjectError e) {    
	cerr << e.getMessage();
    }
    try{
	k = Properties::getIntProperty( "/ExonModel/k" );
    }catch( ProjectError e) {    
	cerr << e.getMessage();
    }
    try{
	patpseudo = Properties::getDoubleProperty( "/ExonModel/patpseudocount" );
    }catch( ProjectError e) { 
	cerr << e.getMessage();
    }
    try{
	exonLenD = Properties::getIntProperty( "/ExonModel/exonlengthD" );
    }catch( ProjectError e) { 
	cerr << e.getMessage();
    }
    try{
	slope_of_bandwidth = Properties::getdoubleProperty( "/ExonModel/slope_of_bandwidth");
    }catch( ProjectError e) { 
	cerr << e.getMessage();
    }
    try{
	minwindowcount = Properties::getIntProperty( "/ExonModel/minwindowcount");
    }catch( ProjectError e) { 
	cerr << e.getMessage();
    }
    try{
	min_exon_length = Properties::getIntProperty( "/ExonModel/minexonlength");
    } catch( ProjectError e) { 
	// optional parameter
    }
    try{
	max_exon_length = Properties::getIntProperty( "/ExonModel/maxexonlength");
    } catch( ProjectError e) { 
	cerr << e.getMessage();
    }
    try{
	minPatSum = Properties::getIntProperty( "/ExonModel/minPatSum");
    } catch( ProjectError e) { 
	cerr << e.getMessage();
    }
    try{
	etorder = Properties::getIntProperty( "/ExonModel/etorder");
    } catch( ProjectError e) { 
	cerr << e.getMessage();
    }
    try{
	etpseudocount = Properties::getIntProperty( "/ExonModel/etpseudocount");
    } catch( ProjectError e) { 
	cerr << e.getMessage();
    }
    try{
	tis_motif_memory = Properties::getIntProperty( "/ExonModel/tis_motif_memory");
    } catch( ProjectError e) { 
	cerr << e.getMessage();
    }
    try{
	tis_motif_radius = Properties::getIntProperty( "/ExonModel/tis_motif_radius");
    } catch( ProjectError e) { 
	cerr << e.getMessage();
    }
	
    trans_init_window = Constant::trans_init_window;

    // reserve space for GC content dependent arrays IF NOT DONE BEFORE
    if (!GCPls) {
      GCPls = new Matrix<vector<Double> >[Constant::decomp_num_steps];
      
      highGCPls = new Matrix<vector<Double> >[Constant::decomp_num_steps];
      mediumGCPls = new Matrix<vector<Double> >[Constant::decomp_num_steps];
      lowGCPls = new Matrix<vector<Double> >[Constant::decomp_num_steps];
    }
    if (!GCtransInitMotif)
      GCtransInitMotif = new Motif[Constant::decomp_num_steps];
    if (!GCetMotif){
      GCetMotif = new Motif**[Constant::decomp_num_steps];
      for (int idx = 0; idx < Constant::decomp_num_steps; idx++){
	GCetMotif[idx] = new Motif*[3];
	GCetMotif[idx][0] = new Motif();
	GCetMotif[idx][1] = new Motif();
	GCetMotif[idx][2] = new Motif();
      }
    }
    if (!GCemiprobs)
      GCemiprobs = new FramedPatMMGroup[Constant::decomp_num_steps];
  
    if (!highGCemiprobs)
      highGCemiprobs = new FramedPatMMGroup[Constant::decomp_num_steps];
    if (!mediumGCemiprobs)
      mediumGCemiprobs = new FramedPatMMGroup[Constant::decomp_num_steps];
    if (!lowGCemiprobs)
      lowGCemiprobs = new FramedPatMMGroup[Constant::decomp_num_steps];
  
    if (!GCinitemiprobs){
      GCinitemiprobs = new vector<Double>*[Constant::decomp_num_steps];
      for (int idx = 0; idx < Constant::decomp_num_steps; idx++)
	GCinitemiprobs[idx] = new vector<Double>[3];
    }

    if (!highGCinitemiprobs){
      highGCinitemiprobs = new vector<Double>*[Constant::decomp_num_steps];
      for (int idx = 0; idx < Constant::decomp_num_steps; idx++)
	highGCinitemiprobs[idx] = new vector<Double>[3];
    }
    if (!mediumGCinitemiprobs){
      mediumGCinitemiprobs = new vector<Double>*[Constant::decomp_num_steps];
      for (int idx = 0; idx < Constant::decomp_num_steps; idx++)
	mediumGCinitemiprobs[idx] = new vector<Double>[3];
    }
    if (!lowGCinitemiprobs){
      lowGCinitemiprobs = new vector<Double>*[Constant::decomp_num_steps];
      for (int idx = 0; idx < Constant::decomp_num_steps; idx++)
	lowGCinitemiprobs[idx] = new vector<Double>[3];
    }

    if (!GCetemiprobs){
      GCetemiprobs = new vector<Double>*[Constant::decomp_num_steps];
      for (int idx = 0; idx < Constant::decomp_num_steps; idx++)
	GCetemiprobs[idx] = new vector<Double>[3];
    }

    if (!highGCetemiprobs){
      highGCetemiprobs = new vector<Double>*[Constant::decomp_num_steps];
      for (int idx = 0; idx < Constant::decomp_num_steps; idx++)
	highGCetemiprobs[idx] = new vector<Double>[3];
    }
    if (!mediumGCetemiprobs){
      mediumGCetemiprobs = new vector<Double>*[Constant::decomp_num_steps];
      for (int idx = 0; idx < Constant::decomp_num_steps; idx++)
	mediumGCetemiprobs[idx] = new vector<Double>[3];
    }
    if (!lowGCetemiprobs){
      lowGCetemiprobs = new vector<Double>*[Constant::decomp_num_steps];
      for (int idx = 0; idx < Constant::decomp_num_steps; idx++)
	lowGCetemiprobs[idx] = new vector<Double>[3];
    }
}

/*
 * ===[ ExonModel::readProbabilities ]====================================
 */
void ExonModel::readProbabilities(int parIndex) {
    if( parIndex == lastParIndex )
	return;
    
    string filename = Constant::fullSpeciesPath() + Properties::getProperty("/ExonModel/infile");
    ifstream istrm; 
    istrm.open(filename.c_str(), ifstream::in);
    if (istrm) {
	int size;
	// Base4Int b4i_e( 0, k+1 );
	Seq2Int s2i_e(k+1);
	int dummyl, dummyi, dummyk;
	Double dbl;

	if (!hasLenDist) {
	    // read length distributions
	    istrm >> goto_line_after( "[LENGTH]");
	    istrm >> comment >> exonLenD;
	    istrm >> comment >> slope_of_bandwidth;
	    istrm >> comment >> minwindowcount;

	    istrm >> comment >> numSingle >> numInitial >> numInternal >> numTerminal;

	    istrm >> comment >> highnumSingle >> highnumInitial >> highnumInternal >> highnumTerminal;
	    istrm >> mediumnumSingle >> mediumnumInitial >> mediumnumInternal >> mediumnumTerminal;
	    istrm >> lownumSingle >> lownumInitial >> lownumInternal >> lownumTerminal;

	    istrm >> comment >> numHugeSingle >> numHugeInitial >> numHugeInternal >> numHugeTerminal;

	    istrm >> comment >> highnumHugeSingle >> highnumHugeInitial >> highnumHugeInternal >> highnumHugeTerminal;
	    istrm >> comment >> mediumnumHugeSingle >> mediumnumHugeInitial >> mediumnumHugeInternal >> mediumnumHugeTerminal;
	    istrm >> comment >> lownumHugeSingle >> lownumHugeInitial >> lownumHugeInternal >> lownumHugeTerminal;

	    istrm >> comment;
	    highlenDistSingle.assign(max_exon_length+1, 0.0);
	    highlenDistInitial.assign(max_exon_length+1, 0.0);
	    highlenDistInternal.assign(max_exon_length+1, 0.0);
	    highlenDistTerminal.assign(max_exon_length+1, 0.0);
	
	    for( int i = 0; i <= exonLenD; i++ ){
		istrm >> dummyi;
		istrm >> dbl;
		highlenDistSingle[i] = dbl / 1000;
		istrm >> dbl; 
		highlenDistInitial[i]= dbl / 1000;
		istrm >> dbl; 
		highlenDistInternal[i] = dbl / 1000;
		istrm >> dbl;
		highlenDistTerminal[i] = dbl / 1000;
	    }
	    // single exon can't be shorter than min_coding_len
	    for ( int i = 0; i < Constant::min_coding_len; i++) 
		highlenDistSingle[i] = 0.0;
	    fillTailsOfLengthDistributionsHigh();

	    istrm >> comment;
	    mediumlenDistSingle.assign(max_exon_length+1, 0.0);
	    mediumlenDistInitial.assign(max_exon_length+1, 0.0);
	    mediumlenDistInternal.assign(max_exon_length+1, 0.0);
	    mediumlenDistTerminal.assign(max_exon_length+1, 0.0);

	    for( int i = 0; i <= exonLenD; i++ ){
		istrm >> dummyi;
		istrm >> dbl;
		mediumlenDistSingle[i] = dbl / 1000;
		istrm >> dbl;
		mediumlenDistInitial[i]= dbl / 1000;
		istrm >> dbl;
		mediumlenDistInternal[i] = dbl / 1000;
		istrm >> dbl;
		mediumlenDistTerminal[i] = dbl / 1000;
	    }
	    // single exon can't be shorter than min_coding_len
	    for ( int i = 0; i < Constant::min_coding_len; i++)
		mediumlenDistSingle[i] = 0.0;
	    fillTailsOfLengthDistributionsMedium();

	    istrm >> comment;
	    lowlenDistSingle.assign(max_exon_length+1, 0.0);
	    lowlenDistInitial.assign(max_exon_length+1, 0.0);
	    lowlenDistInternal.assign(max_exon_length+1, 0.0);
	    lowlenDistTerminal.assign(max_exon_length+1, 0.0);

	    for( int i = 0; i <= exonLenD; i++ ){
		istrm >> dummyi;
		istrm >> dbl;
		lowlenDistSingle[i] = dbl / 1000;
		istrm >> dbl;
		lowlenDistInitial[i]= dbl / 1000;
		istrm >> dbl;
		lowlenDistInternal[i] = dbl / 1000;
		istrm >> dbl;
		lowlenDistTerminal[i] = dbl / 1000;
	    }
	    // single exon can't be shorter than min_coding_len
	    for ( int i = 0; i < Constant::min_coding_len; i++)
		lowlenDistSingle[i] = 0.0;
	    fillTailsOfLengthDistributionsLow();

	    hasLenDist = true;
	}

	    
	/*
	 * begin of content dependent part
	 */

	char zusString[6];
	sprintf(zusString, "[%d]", parIndex);
	istrm >> goto_line_after(zusString);

	/*
	 * content model
	 */

	istrm >> goto_line_after( "[High P_ls]" );
	istrm >> comment;      // dummy k
	highPls.assign(k+1, 3);
	for( int l = 0; l <= k; l++ ){
	    string checkBase;
	    istrm >> comment >> dummyl;
	    int size = POWER4TOTHE(l+1);
	    Seq2Int s2i(l+1);
	    highPls[l][0].resize( size );
	    highPls[l][2] = highPls[l][1] = highPls[l][0];
	    for( int j = 0; j < size; j++ ){
		istrm >> comment;
		int pn = s2i.read(istrm);
		if (pn != j)
		    throw ProjectError("ExonModel::readProbabilities: Error reading file " + filename +
				       " at P_ls, pattern " + s2i.INV(pn));
		istrm >> highPls[l][0][j]
		      >> highPls[l][1][j]
		      >> highPls[l][2][j];
	    }
	}

	istrm >> goto_line_after( "[Medium P_ls]" );
	istrm >> comment;      // dummy k
	mediumPls.assign(k+1, 3);
	for( int l = 0; l <= k; l++ ){
	    string checkBase;
	    istrm >> comment >> dummyl;
	    int size = POWER4TOTHE(l+1);
	    Seq2Int s2i(l+1);
	    mediumPls[l][0].resize( size );
	    mediumPls[l][2] = mediumPls[l][1] = mediumPls[l][0];
	    for( int j = 0; j < size; j++ ){
		istrm >> comment;
		int pn = s2i.read(istrm);
		if (pn != j)
		    throw ProjectError("ExonModel::readProbabilities: Error reading file " + filename +
				       " at P_ls, pattern " + s2i.INV(pn));
		istrm >> mediumPls[l][0][j]
		      >> mediumPls[l][1][j]
		      >> mediumPls[l][2][j];
	    }
	}

	istrm >> goto_line_after( "[Low P_ls]" );
	istrm >> comment;      // dummy k
	lowPls.assign(k+1, 3);
	for( int l = 0; l <= k; l++ ){
	    string checkBase;
	    istrm >> comment >> dummyl;
	    int size = POWER4TOTHE(l+1);
	    Seq2Int s2i(l+1);
	    lowPls[l][0].resize( size );
	    lowPls[l][2] = lowPls[l][1] = lowPls[l][0];
	    for( int j = 0; j < size; j++ ){
		istrm >> comment;
		int pn = s2i.read(istrm);
		if (pn != j)
		    throw ProjectError("ExonModel::readProbabilities: Error reading file " + filename +
				       " at P_ls, pattern " + s2i.INV(pn));
		istrm >> lowPls[l][0][j]
		      >> lowPls[l][1][j]
		      >> lowPls[l][2][j];
	    }
	}

	if (transInitMotif)
	    delete transInitMotif;
	transInitMotif = new Motif();
	istrm >> goto_line_after( "[TRANSINIT]" );
	transInitMotif->read(istrm);

	if (etMotif){
	    delete etMotif[0];
	    delete etMotif[1];
	    delete etMotif[2];
	}

	etMotif = new Motif*[3];
	etMotif[0] = new Motif();
	istrm >> goto_line_after( "[ETMOTIF0]" );
	etMotif[0]->read(istrm);
	etMotif[1] = new Motif();
	istrm >> goto_line_after( "[ETMOTIF1]" );
	etMotif[1]->read(istrm);
	etMotif[2] = new Motif();
	istrm >> goto_line_after( "[ETMOTIF2]" );
	etMotif[2]->read(istrm);

	// EMISSION
	// make the emission probabilities
	for (int f=0; f<3; f++) {
	    highemiprobs.probs[f].resize( highPls[k][f].size() );
	    highemiprobs.order = k;
	}
	// for backward compatibility, check whether EMISSION parameters exist at all, if not compute them from the Psl (old version)
	streampos spos = istrm.tellg();
	istrm >> goto_line_after( "[High EMISSION]" );
	if (!istrm){
	    istrm.clear();
	    istrm.seekg(spos); // go back to where you were
	    for (int f=0; f<3; f++)
		computeEmiFromPat(highPls[k][f], highemiprobs.probs[f], k);
	} else {
	    istrm >> comment >> size;
	    istrm >> comment >> dummyk;
	    istrm >> comment >> dbl;
	    for( int i = 0; i < highemiprobs.probs[0].size(); i++ ){
		istrm >> comment;
		int pn = s2i_e.read(istrm);
		if (pn != i)
		    throw ProjectError("ExonModel::readProbabilities: Error reading file " + filename +
				       " at EMISSION, pattern " + s2i_e.INV(pn));
		istrm >> highemiprobs.probs[0][i]
		      >> highemiprobs.probs[1][i]
		      >> highemiprobs.probs[2][i];
	    }
	}

	// EMISSION
	// make the medium emission probabilities
	for (int f=0; f<3; f++) {
	    mediumemiprobs.probs[f].resize( mediumPls[k][f].size() );
	    mediumemiprobs.order = k;
	}
	// for backward compatibility, check whether EMISSION parameters exist at all, if not compute them from the Psl (old version)
	spos = istrm.tellg();
	istrm >> goto_line_after( "[Medium EMISSION]" );
	if (!istrm){
	    istrm.clear();
	    istrm.seekg(spos); // go back to where you were
	    for (int f=0; f<3; f++)
		computeEmiFromPat(mediumPls[k][f], mediumemiprobs.probs[f], k);
	} else {
	    istrm >> comment >> size;
	    istrm >> comment >> dummyk;
	    istrm >> comment >> dbl;
	    for( int i = 0; i < mediumemiprobs.probs[0].size(); i++ ){
		istrm >> comment;
		int pn = s2i_e.read(istrm);
		if (pn != i)
		    throw ProjectError("ExonModel::readProbabilities: Error reading file " + filename +
				       " at EMISSION, pattern " + s2i_e.INV(pn));
		istrm >> mediumemiprobs.probs[0][i]
		      >> mediumemiprobs.probs[1][i]
		      >> mediumemiprobs.probs[2][i];
	    }
	}

	// EMISSION
	// make the low emission probabilities
	for (int f=0; f<3; f++) {
	    lowemiprobs.probs[f].resize( lowPls[k][f].size() );
	    lowemiprobs.order = k;
	}
	// for backward compatibility, check whether EMISSION parameters exist at all, if not compute them from the Psl (old version)
	spos = istrm.tellg();
	istrm >> goto_line_after( "[Low EMISSION]" );
	if (!istrm){
	    istrm.clear();
	    istrm.seekg(spos); // go back to where you were
	    for (int f=0; f<3; f++)
		computeEmiFromPat(lowPls[k][f], lowemiprobs.probs[f], k);
	} else {
	    istrm >> comment >> size;
	    istrm >> comment >> dummyk;
	    istrm >> comment >> dbl;
	    for( int i = 0; i < lowemiprobs.probs[0].size(); i++ ){
		istrm >> comment;
		int pn = s2i_e.read(istrm);
		if (pn != i)
		    throw ProjectError("ExonModel::readProbabilities: Error reading file " + filename +
				       " at EMISSION, pattern " + s2i_e.INV(pn));
		istrm >> lowemiprobs.probs[0][i]
		      >> lowemiprobs.probs[1][i]
		      >> lowemiprobs.probs[2][i];
	    }
	}
	
	istrm >> goto_line_after( "[High INITEMISSION]" );
	istrm >> comment >> size;
	istrm >> comment >> dummyk;
	if (dummyk != k)
	    throw ProjectError("ExonModel::readProbabilities: Mismatch in order of exon INITEMISSION Markov chain.");
	istrm >> comment >> patpseudo;
	highinitemiprobs[0].resize(size);
	highinitemiprobs[1].resize(size);
	highinitemiprobs[2].resize(size);
	while( istrm >> comment >> ws, istrm && istrm.peek() != '[' ){
	    int pn = s2i_e.read(istrm);
	    istrm >> highinitemiprobs[0][pn]
		  >> highinitemiprobs[1][pn]
		  >> highinitemiprobs[2][pn];
	}

	istrm >> goto_line_after( "[Medium INITEMISSION]" );
	istrm >> comment >> size;
	istrm >> comment >> dummyk;
	if (dummyk != k)
	    throw ProjectError("ExonModel::readProbabilities: Mismatch in order of exon INITEMISSION Markov chain.");
	istrm >> comment >> patpseudo;
	mediuminitemiprobs[0].resize(size);
	mediuminitemiprobs[1].resize(size);
	mediuminitemiprobs[2].resize(size);
	while( istrm >> comment >> ws, istrm && istrm.peek() != '[' ){
	    int pn = s2i_e.read(istrm);
	    istrm >> mediuminitemiprobs[0][pn]
		  >> mediuminitemiprobs[1][pn]
		  >> mediuminitemiprobs[2][pn];
	}

	istrm >> goto_line_after( "[Low INITEMISSION]" );
	istrm >> comment >> size;
	istrm >> comment >> dummyk;
	if (dummyk != k)
	    throw ProjectError("ExonModel::readProbabilities: Mismatch in order of exon INITEMISSION Markov chain.");
	istrm >> comment >> patpseudo;
	lowinitemiprobs[0].resize(size);
	lowinitemiprobs[1].resize(size);
	lowinitemiprobs[2].resize(size);
	while( istrm >> comment >> ws, istrm && istrm.peek() != '[' ){
	    int pn = s2i_e.read(istrm);
	    istrm >> lowinitemiprobs[0][pn]
		  >> lowinitemiprobs[1][pn]
		  >> lowinitemiprobs[2][pn];
	}


	istrm >> goto_line_after( "[High ETEMISSION]" );
	istrm >> comment >> size;
	istrm >> comment >> k;   // TODO: k nicht von hier bestimmen lassen
	istrm >> comment >> patpseudo;
	highetemiprobs[0].resize(size);
	highetemiprobs[1].resize(size);
	highetemiprobs[2].resize(size);
	while( istrm >> comment >> ws, istrm && istrm.peek() != '[' ){
	    int pn = s2i_e.read(istrm);
	    istrm >> highetemiprobs[0][pn]
		  >> highetemiprobs[1][pn]
		  >> highetemiprobs[2][pn];
	}


	istrm >> goto_line_after( "[Medium ETEMISSION]" );
	istrm >> comment >> size;
	istrm >> comment >> k;   // TODO: k nicht von hier bestimmen lassen
	istrm >> comment >> patpseudo;
	mediumetemiprobs[0].resize(size);
	mediumetemiprobs[1].resize(size);
	mediumetemiprobs[2].resize(size);
	while( istrm >> comment >> ws, istrm && istrm.peek() != '[' ){
	    int pn = s2i_e.read(istrm);
	    istrm >> mediumetemiprobs[0][pn]
		  >> mediumetemiprobs[1][pn]
		  >> mediumetemiprobs[2][pn];
	}
	lastParIndex = parIndex;

	istrm >> goto_line_after( "[Low ETEMISSION]" );
	istrm >> comment >> size;
	istrm >> comment >> k;   // TODO: k nicht von hier bestimmen lassen
	istrm >> comment >> patpseudo;
	lowetemiprobs[0].resize(size);
	lowetemiprobs[1].resize(size);
	lowetemiprobs[2].resize(size);
	while( istrm >> comment >> ws, istrm && istrm.peek() != '[' ){
	    int pn = s2i_e.read(istrm);
	    istrm >> lowetemiprobs[0][pn]
		  >> lowetemiprobs[1][pn]
		  >> lowetemiprobs[2][pn];
	}
	lastParIndex = parIndex;
    } else {
	string msg("ExonModel: Couldn't open file ");
	msg += filename;
	throw ProjectError(msg);

    }
}

/*
 * readAllParameters
 * read in species_exon_probs.pbl parameter files and store parameters for all gc content classes
 */
void ExonModel::readAllParameters(){
  string filename = Constant::fullSpeciesPath() + Properties::getProperty("/ExonModel/infile");
  ifstream istrm; 
  istrm.open(filename.c_str(), ifstream::in);
  if (istrm) {
    int size;
    // Base4Int b4i_e( 0, k+1 );
    Seq2Int s2i_e(k+1);
    int dummyl, dummyi, dummyk;
    Double dbl;

    if (!hasLenDist) {
      // read length distributions
      istrm >> goto_line_after( "[LENGTH]");
      istrm >> comment >> exonLenD;
      istrm >> comment >> slope_of_bandwidth;
      istrm >> comment >> minwindowcount;

	    istrm >> comment >> numSingle >> numInitial >> numInternal >> numTerminal;

	    istrm >> comment >> highnumSingle >> highnumInitial >> highnumInternal >> highnumTerminal;
	    istrm >> mediumnumSingle >> mediumnumInitial >> mediumnumInternal >> mediumnumTerminal;
	    istrm >> lownumSingle >> lownumInitial >> lownumInternal >> lownumTerminal;

	    istrm >> comment >> numHugeSingle >> numHugeInitial >> numHugeInternal >> numHugeTerminal;

	    istrm >> comment >> highnumHugeSingle >> highnumHugeInitial >> highnumHugeInternal >> highnumHugeTerminal;
	    istrm >> comment >> mediumnumHugeSingle >> mediumnumHugeInitial >> mediumnumHugeInternal >> mediumnumHugeTerminal;
	    istrm >> comment >> lownumHugeSingle >> lownumHugeInitial >> lownumHugeInternal >> lownumHugeTerminal;

	    istrm >> comment;

      highlenDistSingle.assign(max_exon_length+1, 0.0);
      highlenDistInitial.assign(max_exon_length+1, 0.0);
      highlenDistInternal.assign(max_exon_length+1, 0.0);
      highlenDistTerminal.assign(max_exon_length+1, 0.0);
	
	    for( int i = 0; i <= exonLenD; i++ ){
		istrm >> dummyi;
		istrm >> dbl;
		highlenDistSingle[i] = dbl / 1000;
		istrm >> dbl;
		highlenDistInitial[i]= dbl / 1000;
		istrm >> dbl;
		highlenDistInternal[i] = dbl / 1000;
		istrm >> dbl;
		highlenDistTerminal[i] = dbl / 1000;
	    }
	    // single exon can't be shorter than min_coding_len
	    for ( int i = 0; i < Constant::min_coding_len; i++)
		highlenDistSingle[i] = 0.0;
	    fillTailsOfLengthDistributionsHigh();

	    istrm >> comment;
	    mediumlenDistSingle.assign(max_exon_length+1, 0.0);
	    mediumlenDistInitial.assign(max_exon_length+1, 0.0);
	    mediumlenDistInternal.assign(max_exon_length+1, 0.0);
	    mediumlenDistTerminal.assign(max_exon_length+1, 0.0);

	    for( int i = 0; i <= exonLenD; i++ ){
		istrm >> dummyi;
		istrm >> dbl;
		mediumlenDistSingle[i] = dbl / 1000;
		istrm >> dbl;
		mediumlenDistInitial[i]= dbl / 1000;
		istrm >> dbl;
		mediumlenDistInternal[i] = dbl / 1000;
		istrm >> dbl;
		mediumlenDistTerminal[i] = dbl / 1000;
	    }
	    // single exon can't be shorter than min_coding_len
	    for ( int i = 0; i < Constant::min_coding_len; i++)
		mediumlenDistSingle[i] = 0.0;
	    fillTailsOfLengthDistributionsMedium();

	    istrm >> comment;
	    lowlenDistSingle.assign(max_exon_length+1, 0.0);
	    lowlenDistInitial.assign(max_exon_length+1, 0.0);
	    lowlenDistInternal.assign(max_exon_length+1, 0.0);
	    lowlenDistTerminal.assign(max_exon_length+1, 0.0);

	    for( int i = 0; i <= exonLenD; i++ ){
		istrm >> dummyi;
		istrm >> dbl;
		lowlenDistSingle[i] = dbl / 1000;
		istrm >> dbl;
		lowlenDistInitial[i]= dbl / 1000;
		istrm >> dbl;
		lowlenDistInternal[i] = dbl / 1000;
		istrm >> dbl;
		lowlenDistTerminal[i] = dbl / 1000;
	    }
	    // single exon can't be shorter than min_coding_len
	    for ( int i = 0; i < Constant::min_coding_len; i++)
		lowlenDistSingle[i] = 0.0;
	    fillTailsOfLengthDistributionsLow();

	    hasLenDist = true;

    }


    /*
     * begin of GC content dependent part
     */
    char zusString[6];
    // delete all data structures in case this function has been called before
    /*    if (GCPls)
      delete [] GCPls;
    if (GCtransInitMotif)
	delete [] GCtransInitMotif;
    if (GCetMotif)
      for (int idx = 0; idx < Constant::decomp_num_steps; idx++){
	delete GCetMotif[idx][0];
	delete GCetMotif[idx][1];
	delete GCetMotif[idx][2];
      }
    if (GCemiprobs)
      delete [] GCemiprobs;
    if (GCinitemiprobs)
      for (int idx = 0; idx < Constant::decomp_num_steps; idx++)
      delete [] GCinitemiprobs[idx];*/
      
    // loop over GC content classes and read in all remaining data
    for (int idx = 0; idx < Constant::decomp_num_steps; idx++) {
	  
      sprintf(zusString, "[%d]", idx+1);
      istrm >> goto_line_after(zusString);

      /*
       * content model
       */
	  
      istrm >> goto_line_after( "[High P_ls]" );
      istrm >> comment;      // dummy k 
      highGCPls[idx].assign(k+1, 3);
      for( int l = 0; l <= k; l++ ){
	string checkBase;
	istrm >> comment >> dummyl;
	int size = POWER4TOTHE(l+1);
	Seq2Int s2i(l+1);
	highGCPls[idx][l][0].resize( size );
	highGCPls[idx][l][2] = highGCPls[idx][l][1] = highGCPls[idx][l][0];
	
    
	for( int j = 0; j < size; j++ ){
	  istrm >> comment;
	  int pn = s2i.read(istrm);
	  if (pn != j)
	    throw ProjectError("ExonModel::readProbabilities: Error reading file " + filename +
			       " at P_ls, pattern " + s2i.INV(pn));
	  istrm >> highGCPls[idx][l][0][j]
		>> highGCPls[idx][l][1][j]
		>> highGCPls[idx][l][2][j];
      
	}
      }




        /*
         * content model
         */

        istrm >> goto_line_after( "[Medium P_ls]" );
        istrm >> comment;      // dummy k
        mediumGCPls[idx].assign(k+1, 3);
        for( int l = 0; l <= k; l++ ){
  	string checkBase;
  	istrm >> comment >> dummyl;
  	int size = POWER4TOTHE(l+1);
  	Seq2Int s2i(l+1);
  	mediumGCPls[idx][l][0].resize( size );
  	mediumGCPls[idx][l][2] = mediumGCPls[idx][l][1] = mediumGCPls[idx][l][0];


  	for( int j = 0; j < size; j++ ){
  	  istrm >> comment;
  	  int pn = s2i.read(istrm);
  	  if (pn != j)
  	    throw ProjectError("ExonModel::readProbabilities: Error reading file " + filename +
  			       " at P_ls, pattern " + s2i.INV(pn));
  	  istrm >> mediumGCPls[idx][l][0][j]
  		>> mediumGCPls[idx][l][1][j]
  		>> mediumGCPls[idx][l][2][j];

  	}
        }




          /*
           * content model
           */

          istrm >> goto_line_after( "[Low P_ls]" );
          istrm >> comment;      // dummy k
          lowGCPls[idx].assign(k+1, 3);
          for( int l = 0; l <= k; l++ ){
    	string checkBase;
    	istrm >> comment >> dummyl;
    	int size = POWER4TOTHE(l+1);
    	Seq2Int s2i(l+1);
    	lowGCPls[idx][l][0].resize( size );
    	lowGCPls[idx][l][2] = lowGCPls[idx][l][1] = lowGCPls[idx][l][0];


    	for( int j = 0; j < size; j++ ){
    	  istrm >> comment;
    	  int pn = s2i.read(istrm);
    	  if (pn != j)
    	    throw ProjectError("ExonModel::readProbabilities: Error reading file " + filename +
    			       " at P_ls, pattern " + s2i.INV(pn));
    	  istrm >> lowGCPls[idx][l][0][j]
    		>> lowGCPls[idx][l][1][j]
    		>> lowGCPls[idx][l][2][j];

    	}
          }


      istrm >> goto_line_after( "[TRANSINIT]" );
      GCtransInitMotif[idx].read(istrm);
      istrm >> goto_line_after( "[ETMOTIF0]" );
      GCetMotif[idx][0]->read(istrm);
      istrm >> goto_line_after( "[ETMOTIF1]" );
      GCetMotif[idx][1]->read(istrm);
      istrm >> goto_line_after( "[ETMOTIF2]" );
      GCetMotif[idx][2]->read(istrm);

      // EMISSION
      // make the emission probabilities
      
      highGCemiprobs[idx].setName("high exon emiprob gc" + (idx+1));
      for (int f=0; f<3; f++) {
	highGCemiprobs[idx].probs[f].resize( highGCPls[idx][k][f].size() );
	highGCemiprobs[idx].order = k;
      }  

      mediumGCemiprobs[idx].setName("medium exon emiprob gc" + (idx+1));
      for (int f=0; f<3; f++) {
	mediumGCemiprobs[idx].probs[f].resize( mediumGCPls[idx][k][f].size() );
	mediumGCemiprobs[idx].order = k;
      }    

      lowGCemiprobs[idx].setName("low exon emiprob gc" + (idx+1));
      for (int f=0; f<3; f++) {
	lowGCemiprobs[idx].probs[f].resize( lowGCPls[idx][k][f].size() );
	lowGCemiprobs[idx].order = k;
      }      
      // for backward compatibility, check whether EMISSION parameters exist at all, if not compute them from the Psl (old version)
      streampos spos = istrm.tellg();
      istrm >> goto_line_after( "[High EMISSION]" );
      if (!istrm){
	istrm.clear();
	istrm.seekg(spos); // go back to where you were
	for (int f=0; f<3; f++) {
	  computeEmiFromPat(highGCPls[idx][k][f], highGCemiprobs[idx].probs[f], k);
      }
    } else {
	istrm >> comment >> size;
	istrm >> comment >> dummyk;
	istrm >> comment >> dbl;
	for( int i = 0; i < highGCemiprobs[idx].probs[0].size(); i++ ){
	  istrm >> comment;
	  int pn = s2i_e.read(istrm);
	  if (pn != i)
	    throw ProjectError("ExonModel::readProbabilities: Error reading file " + filename +
			       " at EMISSION, pattern " + s2i_e.INV(pn));
	  istrm >> highGCemiprobs[idx].probs[0][i]
		>> highGCemiprobs[idx].probs[1][i]
		>> highGCemiprobs[idx].probs[2][i];
	}
      }

      spos = istrm.tellg();
      istrm >> goto_line_after( "[Medium EMISSION]" );
      if (!istrm){
	istrm.clear();
	istrm.seekg(spos); // go back to where you were
	for (int f=0; f<3; f++) {
	  computeEmiFromPat(mediumGCPls[idx][k][f], mediumGCemiprobs[idx].probs[f], k);
      }
    } else {
	istrm >> comment >> size;
	istrm >> comment >> dummyk;
	istrm >> comment >> dbl;
	for( int i = 0; i < mediumGCemiprobs[idx].probs[0].size(); i++ ){
	  istrm >> comment;
	  int pn = s2i_e.read(istrm);
	  if (pn != i)
	    throw ProjectError("ExonModel::readProbabilities: Error reading file " + filename +
			       " at EMISSION, pattern " + s2i_e.INV(pn));
	  istrm >> mediumGCemiprobs[idx].probs[0][i]
		>> mediumGCemiprobs[idx].probs[1][i]
		>> mediumGCemiprobs[idx].probs[2][i];
	}
      }


      spos = istrm.tellg();
      istrm >> goto_line_after( "[Low EMISSION]" );
      if (!istrm){
	istrm.clear();
	istrm.seekg(spos); // go back to where you were
	for (int f=0; f<3; f++) {
	  computeEmiFromPat(lowGCPls[idx][k][f], lowGCemiprobs[idx].probs[f], k);
      }
    } else {
	istrm >> comment >> size;
	istrm >> comment >> dummyk;
	istrm >> comment >> dbl;
	for( int i = 0; i < lowGCemiprobs[idx].probs[0].size(); i++ ){
	  istrm >> comment;
	  int pn = s2i_e.read(istrm);
	  if (pn != i)
	    throw ProjectError("ExonModel::readProbabilities: Error reading file " + filename +
			       " at EMISSION, pattern " + s2i_e.INV(pn));
	  istrm >> lowGCemiprobs[idx].probs[0][i]
		>> lowGCemiprobs[idx].probs[1][i]
		>> lowGCemiprobs[idx].probs[2][i];
	}
      }

      istrm >> goto_line_after( "[High INITEMISSION]" );
      istrm >> comment >> size;
      istrm >> comment >> dummyk;
      if (dummyk != k)
	throw ProjectError("ExonModel::readProbabilities: Mismatch in order of exon INITEMISSION Markov chain.");
      istrm >> comment >> patpseudo;
      highGCinitemiprobs[idx][0].resize(size);
      highGCinitemiprobs[idx][1].resize(size);
      highGCinitemiprobs[idx][2].resize(size);
      while( istrm >> comment >> ws, istrm && istrm.peek() != '[' ){
	int pn = s2i_e.read(istrm);
	istrm >> highGCinitemiprobs[idx][0][pn]
	      >> highGCinitemiprobs[idx][1][pn]
	      >> highGCinitemiprobs[idx][2][pn];
      }

      istrm >> goto_line_after( "[Medium INITEMISSION]" );
      istrm >> comment >> size;
      istrm >> comment >> dummyk;
      if (dummyk != k)
	throw ProjectError("ExonModel::readProbabilities: Mismatch in order of exon INITEMISSION Markov chain.");
      istrm >> comment >> patpseudo;
      mediumGCinitemiprobs[idx][0].resize(size);
      mediumGCinitemiprobs[idx][1].resize(size);
      mediumGCinitemiprobs[idx][2].resize(size);
      while( istrm >> comment >> ws, istrm && istrm.peek() != '[' ){
	int pn = s2i_e.read(istrm);
	istrm >> mediumGCinitemiprobs[idx][0][pn]
	      >> mediumGCinitemiprobs[idx][1][pn]
	      >> mediumGCinitemiprobs[idx][2][pn];
      }

      istrm >> goto_line_after( "[Low INITEMISSION]" );
      istrm >> comment >> size;
      istrm >> comment >> dummyk;
      if (dummyk != k)
	throw ProjectError("ExonModel::readProbabilities: Mismatch in order of exon INITEMISSION Markov chain.");
      istrm >> comment >> patpseudo;
      lowGCinitemiprobs[idx][0].resize(size);
      lowGCinitemiprobs[idx][1].resize(size);
      lowGCinitemiprobs[idx][2].resize(size);
      while( istrm >> comment >> ws, istrm && istrm.peek() != '[' ){
	int pn = s2i_e.read(istrm);
	istrm >> lowGCinitemiprobs[idx][0][pn]
	      >> lowGCinitemiprobs[idx][1][pn]
	      >> lowGCinitemiprobs[idx][2][pn];
      }

      istrm >> goto_line_after( "[High ETEMISSION]" );
      istrm >> comment >> size;
      istrm >> comment >> k;
      istrm >> comment >> patpseudo;
      highGCetemiprobs[idx][0].resize(size);
      highGCetemiprobs[idx][1].resize(size);
      highGCetemiprobs[idx][2].resize(size);
      while( istrm >> comment >> ws, istrm && istrm.peek() != '[' ){
	int pn = s2i_e.read(istrm);
	istrm >> highGCetemiprobs[idx][0][pn]
	      >> highGCetemiprobs[idx][1][pn]
	      >> highGCetemiprobs[idx][2][pn];
      }

      istrm >> goto_line_after( "[Medium ETEMISSION]" );
      istrm >> comment >> size;
      istrm >> comment >> k;
      istrm >> comment >> patpseudo;
      mediumGCetemiprobs[idx][0].resize(size);
      mediumGCetemiprobs[idx][1].resize(size);
      mediumGCetemiprobs[idx][2].resize(size);
      while( istrm >> comment >> ws, istrm && istrm.peek() != '[' ){
	int pn = s2i_e.read(istrm);
	istrm >> mediumGCetemiprobs[idx][0][pn]
	      >> mediumGCetemiprobs[idx][1][pn]
	      >> mediumGCetemiprobs[idx][2][pn];
      }

      istrm >> goto_line_after( "[Low ETEMISSION]" );
      istrm >> comment >> size;
      istrm >> comment >> k;
      istrm >> comment >> patpseudo;
      lowGCetemiprobs[idx][0].resize(size);
      lowGCetemiprobs[idx][1].resize(size);
      lowGCetemiprobs[idx][2].resize(size);
      while( istrm >> comment >> ws, istrm && istrm.peek() != '[' ){
	int pn = s2i_e.read(istrm);
	istrm >> lowGCetemiprobs[idx][0][pn]
	      >> lowGCetemiprobs[idx][1][pn]
	      >> lowGCetemiprobs[idx][2][pn];
      }
    }
    istrm.close();
  } else {
    string msg("ExonModel: Couldn't open file ");
    msg += filename;
    throw ProjectError(msg);

  }
}

/*
 * ===[ ExonModel::fillTailsOfLengthDistributions]====================================
 * Define the tail of the length distributions: from exonLenD up to max_exon_length
 * Make it so that there is no jump at exonLenD and so that the probability of huge
 * exons (>= d) is correct.
 * Probabilities for lengths greater than max_exon_length are 0.
 */
void ExonModel::fillTailsOfLengthDistributions( ){
    Double a,p; 
    int k;
    
    a = lenDistSingle[exonLenD];
    p = Double(1.0) - a/(Double(numHugeSingle+1)/Double(numSingle+1));
    for (k = exonLenD+1; k <= max_exon_length; k++){
	lenDistSingle[k] = p * lenDistSingle[k-1];
    } 

    a = lenDistInitial[exonLenD];
    p = Double(1.0) - a/((Double(numHugeInitial)+1)/(numInitial+1));
    for (k = exonLenD+1; k <= max_exon_length; k++){
	lenDistInitial[k] = p * lenDistInitial[k-1];
    } 

    a = lenDistInternal[exonLenD];
    p = Double(1.0) - a/((Double(numHugeInternal)+1)/(numInternal+1));
    for (k = exonLenD+1; k <= max_exon_length; k++){
	lenDistInternal[k] = p * lenDistInternal[k-1];
    } 

    a = lenDistTerminal[exonLenD];
    p = Double(1.0) - a/((Double(numHugeTerminal)+1)/(numTerminal+1));
    for (k = exonLenD+1; k <= max_exon_length; k++){
	lenDistTerminal[k] = p * lenDistTerminal[k-1];
    }
}

void ExonModel::fillTailsOfLengthDistributionsHigh( ){
    Double a,p;
    int k;

    a = highlenDistSingle[exonLenD];
    p = Double(1.0) - a/(Double(highnumHugeSingle+1)/Double(highnumSingle+1));
    for (k = exonLenD+1; k <= max_exon_length; k++){
	highlenDistSingle[k] = p * highlenDistSingle[k-1];
    }

    a = highlenDistInitial[exonLenD];
    p = Double(1.0) - a/((Double(highnumHugeInitial)+1)/(highnumInitial+1));
    for (k = exonLenD+1; k <= max_exon_length; k++){
	highlenDistInitial[k] = p * highlenDistInitial[k-1];
    }

    a = highlenDistInternal[exonLenD];
    p = Double(1.0) - a/((Double(highnumHugeInternal)+1)/(highnumInternal+1));
    for (k = exonLenD+1; k <= max_exon_length; k++){
	highlenDistInternal[k] = p * highlenDistInternal[k-1];
    }

    a = highlenDistTerminal[exonLenD];
    p = Double(1.0) - a/((Double(highnumHugeTerminal)+1)/(highnumTerminal+1));
    for (k = exonLenD+1; k <= max_exon_length; k++){
	highlenDistTerminal[k] = p * highlenDistTerminal[k-1];
    }
}

void ExonModel::fillTailsOfLengthDistributionsMedium( ){
    Double a,p;
    int k;

    a = mediumlenDistSingle[exonLenD];
    p = Double(1.0) - a/(Double(mediumnumHugeSingle+1)/Double(mediumnumSingle+1));
    for (k = exonLenD+1; k <= max_exon_length; k++){
	mediumlenDistSingle[k] = p * mediumlenDistSingle[k-1];
    }

    a = mediumlenDistInitial[exonLenD];
    p = Double(1.0) - a/((Double(mediumnumHugeInitial)+1)/(mediumnumInitial+1));
    for (k = exonLenD+1; k <= max_exon_length; k++){
	mediumlenDistInitial[k] = p * mediumlenDistInitial[k-1];
    }

    a = mediumlenDistInternal[exonLenD];
    p = Double(1.0) - a/((Double(mediumnumHugeInternal)+1)/(mediumnumInternal+1));
    for (k = exonLenD+1; k <= max_exon_length; k++){
	mediumlenDistInternal[k] = p * mediumlenDistInternal[k-1];
    }

    a = mediumlenDistTerminal[exonLenD];
    p = Double(1.0) - a/((Double(mediumnumHugeTerminal)+1)/(mediumnumTerminal+1));
    for (k = exonLenD+1; k <= max_exon_length; k++){
	mediumlenDistTerminal[k] = p * mediumlenDistTerminal[k-1];
    }
}

void ExonModel::fillTailsOfLengthDistributionsLow( ){
    Double a,p;
    int k;

    a = lowlenDistSingle[exonLenD];
    p = Double(1.0) - a/(Double(lownumHugeSingle+1)/Double(lownumSingle+1));
    for (k = exonLenD+1; k <= max_exon_length; k++){
	lowlenDistSingle[k] = p * lowlenDistSingle[k-1];
    }

    a = lowlenDistInitial[exonLenD];
    p = Double(1.0) - a/((Double(lownumHugeInitial)+1)/(lownumInitial+1));
    for (k = exonLenD+1; k <= max_exon_length; k++){
	lowlenDistInitial[k] = p * lowlenDistInitial[k-1];
    }

    a = lowlenDistInternal[exonLenD];
    p = Double(1.0) - a/((Double(lownumHugeInternal)+1)/(lownumInternal+1));
    for (k = exonLenD+1; k <= max_exon_length; k++){
	lowlenDistInternal[k] = p * lowlenDistInternal[k-1];
    }

    a = lowlenDistTerminal[exonLenD];
    p = Double(1.0) - a/((Double(lownumHugeTerminal)+1)/(lownumTerminal+1));
    for (k = exonLenD+1; k <= max_exon_length; k++){
	lowlenDistTerminal[k] = p * lowlenDistTerminal[k-1];
    }
}


/*
 * ===[ ExonModel::initAlgorithms ]=======================================
 */
void ExonModel::initAlgorithms(Matrix<Double>& trans, int cur) {
  if (!initAlgorithmsCalled) {
    // stuff that needs to be called once for all exon states
    // set these parameters to the one of the GC content index gcIdx    
    Pls = GCPls[gcIdx];

    highPls = highGCPls[gcIdx];
    mediumPls = mediumGCPls[gcIdx];
    lowPls = lowGCPls[gcIdx];

    emiprobs = GCemiprobs[gcIdx];

    highemiprobs = highGCemiprobs[gcIdx];
    mediumemiprobs = mediumGCemiprobs[gcIdx];
    lowemiprobs = lowGCemiprobs[gcIdx];

    initemiprobs[0] =  GCinitemiprobs[gcIdx][0];

    highinitemiprobs[0] =  highGCinitemiprobs[gcIdx][0];
    mediuminitemiprobs[0] =  mediumGCinitemiprobs[gcIdx][0];
    lowinitemiprobs[0] =  lowGCinitemiprobs[gcIdx][0];

    initemiprobs[1] =  GCinitemiprobs[gcIdx][1];
    
    highinitemiprobs[1] =  highGCinitemiprobs[gcIdx][1];
    mediuminitemiprobs[1] =  mediumGCinitemiprobs[gcIdx][1];
    lowinitemiprobs[1] =  lowGCinitemiprobs[gcIdx][1];
    
   
    initemiprobs[2] =  GCinitemiprobs[gcIdx][2];
    
    highinitemiprobs[2] =  highGCinitemiprobs[gcIdx][2];
    mediuminitemiprobs[2] =  mediumGCinitemiprobs[gcIdx][2];
    lowinitemiprobs[2] =  lowGCinitemiprobs[gcIdx][2];
    
    etemiprobs[0] =  GCetemiprobs[gcIdx][0];
    
    highetemiprobs[0] =  highGCetemiprobs[gcIdx][0];
    mediumetemiprobs[0] =  mediumGCetemiprobs[gcIdx][0];
    lowetemiprobs[0] =  lowGCetemiprobs[gcIdx][0];
    
    etemiprobs[1] =  GCetemiprobs[gcIdx][1];
    
    highetemiprobs[1] =  highGCetemiprobs[gcIdx][1];
    mediumetemiprobs[1] =  mediumGCetemiprobs[gcIdx][1];
    lowetemiprobs[1] =  lowGCetemiprobs[gcIdx][1];
    
    etemiprobs[2] =  GCetemiprobs[gcIdx][2];
    
    highetemiprobs[2] =  highGCetemiprobs[gcIdx][2];
    mediumetemiprobs[2] =  mediumGCetemiprobs[gcIdx][2];
    lowetemiprobs[2] =  lowGCetemiprobs[gcIdx][2];
    
    transInitMotif = &GCtransInitMotif[gcIdx];
    etMotif = GCetMotif[gcIdx];
  }
  initAlgorithmsCalled = true;
}


/*
 * ===[ ExonModel::viterbiForwardAndSampling ]=====================================
 *
 * Find the predecessor state "predState" and the ending position "predBase" of 
 * this predecessor state in a most probable path ending in state "state" at position "base". 
 *
 */
void ExonModel::viterbiForwardAndSampling(ViterbiMatrixType& viterbi, // matrix of viterbi variables
					  ViterbiMatrixType& forward, // matrix of forward variables
					  int state,
					  int base,                   // end of the exon state
					  AlgorithmVariant algovar,   // doViterbiOnly, doViterbiAndForward, doSampling or doBacktracking
					  OptionListItem& oli) const {
/* 
 *  an example:
 *  pred.State |ATG               | <= (k+1)     | *** ... ***   | STP        |
 *             |start codon       | pattern prob | emission prob | stop codon |
 *             |                  |              | >= 0          |            |
 *                                 <---- inner sequence  ------->
 *  --------------------------------------------------------------------------
 *  predProb   |beginPartProb     | startSeqProb | restSeqProb   | endPartProb|
 *              <--------------------- lenPartProb -------------------------->
 * 
 */ 

    Double fwdsum = 0;
    Double maxProb = 0;
    OptionsList *optionslist = NULL;
    PP::ExonScorer* exonScorer = 0;
    if (algovar==doSampling)
	optionslist = new OptionsList();
    bool checkSubstates = (profileModel != 0); // && isOnFStrand(etype) && etype != singleG

    /*
     * endPartProb
     * probability of the last - fixed length - part of the exon:
     * translation termination if the exon is a terminal exon
     */
    Double endPartProb = endPartEmiProb(base);
    ViterbiSubmapType substates(state);
    
    /*
     * right is the rightmost position of the inner sequence part
     * just in front of the endPart
     * endOfBioExon is the position most downstream (wrt forward strand) of the biological exon
     */
    
    int endOfBioExon = base + baseOffset; 
    int right = endOfBioExon - innerPartEndOffset;
    
    // if (etype == singleG || etype == terminal) {
    // 	right = base - STOPCODON_LEN;
    // 	endOfBioExon = base;
    // } else if(etype == rsingleG || etype == rinitial) {
    // 	right = base - trans_init_window - STARTCODON_LEN;
    // 	endOfBioExon = base - trans_init_window;
    // } else if(isOnFStrand(etype)){
    // 	right = base;
    // 	endOfBioExon = base + Constant::dss_start;
    // } else {
    // 	right = base;
    // 	endOfBioExon = base + Constant::ass_end;
    // }
    
    /*
     * If it is impossible that a (complete) exon ends here return immediately.
     */
    
    if (endPartProb <= 0 || right < 0) {
	// viterbi[base].erase(state);
	// if (needForwardTable(algovar))
	//     forward[base].erase(state);
	return;
    }

    /*
     * Reading frame of the position "right"
     * win=0,1,2 is the reading frame of __ endOfBioExon+1 __
     * on the reverse strand, count backwards 
     */
    int frameOfRight = isOnFStrand(etype) ? 
	mod3(win - (endOfBioExon + 1) + right) : 
	mod3(win +  endOfBioExon + 1  - right);
    
    // Need this for protein pattern evaluation
    // int endOfLastCodon = endOfBioExon - (isOnFStrand(etype) ? win : 2 - win);
    // if (endOfLastCodon >= dnalen)
    // 	checkSubstates=false;
    if (endOfBioExon >= dnalen)
	checkSubstates = false;

     /*
     * Determine the leftmost possible boundary of the inner sequence part.  
     * That is, determine the longest possible open reading frame.
     */
    int endOfNonStopcodon = (etype == terminalHigh || etype == terminalMedium || etype == terminalLow || etype == singleGHigh || etype == singleGMedium || etype == singleGLow) ?
	endOfBioExon - STOPCODON_LEN : endOfBioExon;
    if (endOfNonStopcodon > dnalen-1) // to be safe in those rare cases
	endOfNonStopcodon = dnalen-1;

    int frameOfEndOfNonStopcodon = isOnFStrand(etype) ?
	mod3(win - 1 - endOfBioExon + endOfNonStopcodon) : 
	mod3(win + 1 + endOfBioExon - endOfNonStopcodon);

   // switch(etype) { 
    // 	case singleG: case initial0: case initial1: case initial2:
    // 	    left = ORFleft + STARTCODON_LEN; break;
    // 	case terminal: case internal0: case internal1: case internal2:
    // 	  if (ORFleft > 0)
    // 	    left = ORFleft + Constant::ass_end; 
    // 	  else
    // 	    left = 0; 
    // 	  break;
    //    	case rinitial: case rinternal0: case rinternal1: case rinternal2:
    // 	  if (ORFleft > 0)
    // 	    left = ORFleft + Constant::dss_start;
    // 	  else
    // 	    left = 0;
    // 	  break;
    // 	default: // i.e., when there must be a reverse stop codon at the left end ...
    // 	    left = ORFleft;
    // }
    
    int ORFleft = orf->leftmostExonBegin(frameOfEndOfNonStopcodon, endOfNonStopcodon, isOnFStrand(etype));

    /*
     * get the extrinsic exonpart information about parts falling in this range
     * and with fitting reading frame
     */
    int seqRelFrame = isOnFStrand(etype)? mod3(right - frameOfRight) : mod3(right + frameOfRight);
    Feature* extrinsicexons = seqFeatColl->getExonListOvlpingRange(ORFleft-1, endOfBioExon, bothstrands, seqRelFrame);

    if (checkSubstates) {
	if (algovar == doBacktracking) {
	    // in backtracking, state is a FullStateId if etype!=terminal
	    SubstateId substate;
	    getStatePair(state, state, substate); // extract the main state index
#ifdef DEBUG
#ifndef DEBUG_BACK
	    if (!substate.undef())
		cerr << "substate bonus at base " << base << ": "
		     << viterbi[base].get(state, substate) / viterbi[base][state]
		     << endl;
#endif
#endif		    
	    // if it's a terminal or rinitial exon, we backtrack from igenic which is
	    // unaware of substates, so the given state won't refer to a substate; 
	    // in order to find out, we check viterbi: if the probability comes from the
	    // profile, there will be substates, otherwise not
	    if (isLastExon(etype))  {
#ifdef DEBUG
		if (!substate.undef())
		    throw ProjectError("viterbiForwardAndSampling called with"
				       " bad state/substate combination"
			               " (last exon but substate=" + itoa(substate) + ").");
#endif
		checkSubstates = viterbi[base].hasSubstates(state);
	    } else {
		checkSubstates = !substate.undef();
	    }
	    if (checkSubstates) 
		exonScorer = new PP::SingleTargetExonScorer(*profileModel, etype, substate, 
							    endOfBioExon, ORFleft);
	} else {
	    // profileModel->advanceHSColl(endOfLastCodon);
	    if (isLastExon(etype))
		exonScorer = new PP::SingleTargetExonScorer(*profileModel, etype,
							    endOfBioExon);
	    else {
		substates.activate();
		exonScorer = new PP::MultiTargetExonScorer(*profileModel, substates, etype,
							   endOfBioExon);
	    }
	}
    }
#ifdef DEBUG
    if (algovar==doBacktracking)
	oli.base = oli.state = -1;
#endif

    // determine bounds for exon start
    int startMax = endOfBioExon + innerPartOffset - min_exon_length + 1;
    int startMin;
    if (isRTerminalExon(etype) || etype == rsingleGHigh || etype == rsingleGMedium || etype == rsingleGLow) 
	startMin = startMax = ORFleft+2; // single gene start candidate position
    else {
	if (ORFleft <= 0)
	    startMin = 0; // left truncated exon
	else
	    // innerPartOffset is beginOfStart - beginOfBioExon
	    startMin = ORFleft + innerPartOffset; 
	if (startMax > base + beginPartLen) // ensure that base>endOfPred
	    startMax = base + beginPartLen;
    }
     
    /*
     * loop over the length of the inner sequence
     */
    seqProb(-1,0,0);      // initialize the static variables
    for (int beginOfStart = startMax; beginOfStart >= startMin; beginOfStart--) {
	// main work done in this funcion call
	int endOfPred = beginOfStart - beginPartLen -1;
	Double notEndPartProb = notEndPartEmiProb(beginOfStart, right, frameOfRight, extrinsicexons);
	if (notEndPartProb <= 0.0 || endOfPred >= dnalen) 
	    continue;
	// one of two cases where endOfPred < 0
	// - initial or single exon starts right at pos 1
	// - left truncated internal or terminal exon
	const ViterbiColumnType& predForw =  forward[endOfPred >= 0 ? endOfPred : 0];
	ViterbiColumnType& predVit = viterbi[endOfPred >= 0 ? endOfPred : 0];
	// compute the maximum over the predecessor states probs times transition probability
	int predReadingFrame;
	int beginOfBioExon  = beginOfStart - innerPartOffset;
	int exonLength = endOfBioExon - beginOfBioExon +1;
	if (checkSubstates)
	    exonScorer->newFirstCodon(beginOfBioExon);
	
	vector<Ancestor>::const_iterator it;
	for( it = ancestor.begin(); it != ancestor.end(); ++it ){
	    int predState = it->pos;
	    if (algovar == doSampling) {
		if (predForw[predState] == 0)
		    continue;
	    } else
		if (predVit.get(predState)==0)
		    continue;
	    StateType predStateType = (*stateMap)[predState];
	    Double transEmiProb = it->val * endPartProb * notEndPartProb;
	    predReadingFrame = stateReadingFrames[predStateType];
	    if (etype == singleGHigh || etype == singleGMedium || etype == singleGLow || 
	    	etype == rsingleGHigh || etype == rsingleGMedium || etype == rsingleGLow ||
	    	etype == rterminal0High || etype == rterminal0Medium || etype == rterminal0Low ||
	    	etype == rterminal1High || etype == rterminal1Medium || etype == rterminal1Low || 
		etype == rterminal2High || etype == rterminal2Medium || etype == rterminal2Low ||
		etype == initial0High || etype == initial0Medium || etype == initial0Low ||
		etype == initial1High || etype == initial1Medium || etype == initial1Low ||
		etype == initial2High || etype == initial2Medium || etype == initial2Low ||
		// has been checked in notEndPartProb
		win == mod3((isOnFStrand(etype))? predReadingFrame + exonLength: predReadingFrame - exonLength)) 
	    {
		// transition and exon length fit with respect to reading frame
		if (needForwardTable(algovar)) {
		    Double fwdsummand = predForw[predState] * transEmiProb;
		    if (algovar == doSampling) {
			if (fwdsummand > 0) 
			    optionslist->add(predState, endOfPred, fwdsummand);
			continue; // next for
		    } else 
			fwdsum += fwdsummand;
		}
		if (checkSubstates) {
		    if (exonScorer->validSize())
			try {
			    exonScorer->score(predVit, predState, transEmiProb);
			    if (algovar==doBacktracking && exonScorer->hasNewMax()) {
				oli.base = endOfPred;
				oli.state = getFullStateId(predState, exonScorer->getPredSubstate());
			    }
			} catch(NoSubmapFoundError e) {}
		} 
		Double predProb = predVit.get(predState) * transEmiProb;
		if (predProb > maxProb) {
		    maxProb = predProb;
		    if (algovar==doBacktracking && !checkSubstates) {
			oli.base = endOfPred;
			oli.state = predState;
		    }
		}
	    }
	} // end of loop over predecessor state
    } // end of loop over the exon length
    
    switch (algovar) {
	case doSampling:
	    optionslist->prepareSampling();
	    try {
		oli = optionslist->sample();
	    } catch (ProjectError e) {
		cerr << "Sampling error in exon model. state=" << state << " base=" << base << endl;
		throw e;
	    }
	    delete optionslist;
	    break;
	case doViterbiAndForward:
	    if (fwdsum > 0)
		forward[base][state] = fwdsum;
	case doViterbiOnly:
	    if (checkSubstates) {
		exonScorer->postProcessing(maxProb);
		if (isLastExon(etype)) {
		    substates.activate();
		    exonScorer->exportSubstates(substates);
		    substates.useMaximum(maxProb);
		}
		if (!substates.empty()) {
		    if (!isFirstExon(etype))
			substates.registerPredecessors();
		    viterbi[base].addSubstates(substates);
		}
	    }
	    if (maxProb > 0)
		viterbi[base][state] = maxProb;
	    break;
	case doBacktracking:
	    if (checkSubstates) {
		exonScorer->addMatches(oli.base + beginPartLen + 1 - innerPartOffset, endOfBioExon);
#ifdef DEBUG
		static_cast<PP::SingleTargetExonScorer*>(exonScorer)->
		    checkResults(oli, viterbi, state, base, endPartProb, 
				 oli.state==-1 ? 0 : notEndPartEmiProb(oli.base + beginPartLen + 1, right, 
								       frameOfRight, extrinsicexons));
#endif			
	    }
	default:
	    break;
    }
    delete exonScorer;
} // ExonModel::viterbiForwardAndSampling


/*
 * ===[ ExonModel::endPartEmiProb ]=====================================
 *
 * Compute the probability that the end part ends at "end" in the test sequence.
 * Return 0 if it is impossible.   
 */

Double ExonModel::endPartEmiProb(int end) const {
    static Double endPartProb = 0.0;
    Feature *feature;
    Double extrinsicEmiQuot = 1.0;
    switch( etype ){
        case singleGHigh: case singleGMedium: case singleGLow: case terminalHigh: case terminalMedium: case terminalLow:
  	{
	    int stppos = end - STOPCODON_LEN + 1;
	    if( stppos < 0 || stppos > dnalen-3 || !GeneticCode::isStopcodon(sequence + stppos) )
	      endPartProb = 0.0;
	    else { 
		// assign probabilities to the stop codons
		// human: taa 28%, tga 48%, tag 24%
		if (strncmp(sequence + stppos, "taa", 3)==0)
		    endPartProb = Constant::ochreprob;
		else if (strncmp(sequence + stppos, "tag", 3)==0)
		    endPartProb = Constant::amberprob;
		else if (strncmp(sequence + stppos, "tga", 3)==0)
		    endPartProb = Constant::opalprob;
		else 
		    throw ProjectError("ExonModel::endPartEmiProb: internal error, unknown stop codon");
		// check if we have extrinsic information about a stop codon
		feature = seqFeatColl->getFeatureListOvlpingRange(stopF, end-2, end , plusstrand);
		if (feature) {
		  while (feature) {
		    if (feature->start <= end-2 && feature->end >= end)
		      extrinsicEmiQuot *= feature->distance_faded_bonus(end-1);
		    feature = feature->next;
		  }
		} else if (seqFeatColl->collection->hasHintsFile)
		  extrinsicEmiQuot = seqFeatColl->collection->malus(stopF);
		/*
		feature = seqFeatColl->getFeatureAt(stopF, end, plusstrand); 
		if (feature)
		    extrinsicEmiQuot = feature->bonus;
		else if (seqFeatColl->collection->hasHintsFile) {
		    extrinsicEmiQuot = seqFeatColl->collection->malus(stopF);
		    }*/
	    }
	}
	break;
	case rsingleGHigh: case rsingleGMedium: case rsingleGLow: case rinitialHigh: case rinitialMedium: case rinitialLow:
	{
	    int startpos = end - trans_init_window - STARTCODON_LEN + 1;
	    if (startpos >= 0 && onRStart(sequence + startpos)) {
		if (startpos + STARTCODON_LEN + trans_init_window - 1 + tis_motif_memory < dnalen)
		    endPartProb = transInitMotif->seqProb(sequence + startpos + STARTCODON_LEN, true, true);  
		else
		    endPartProb = pow(0.25, (double)(dnalen-(startpos + STARTCODON_LEN)));
	    } else 
		endPartProb = 0.0;
	    // check if we have extrinsic information about a reverse start codon
	    feature = seqFeatColl->getFeatureListOvlpingRange(startF, startpos, startpos + STARTCODON_LEN - 1 , minusstrand); 
	    if (feature) {
	      while (feature) {
		if (feature->start <= startpos && feature->end>= startpos + STARTCODON_LEN - 1)
		  extrinsicEmiQuot *= feature->distance_faded_bonus(startpos + 1);
		feature = feature->next;
	      }
	    } else if (seqFeatColl->collection->hasHintsFile)
	      extrinsicEmiQuot = seqFeatColl->collection->malus(startF);
	      
	    //feature = seqFeatColl->getFeatureAt(startF, startpos + STARTCODON_LEN - 1 , minusstrand); 
	    //if (feature)
	    //	extrinsicEmiQuot = feature->bonus;
	    //else if (seqFeatColl->collection->hasHintsFile)
	    //	extrinsicEmiQuot = seqFeatColl->collection->malus(startF);
	}
	break;
        case initial0High: case initial0Medium: case initial0Low: 
        case initial1High: case initial1Medium: case initial1Low: 
        case initial2High: case initial2Medium: case initial2Low:
        case internal0High: case internal0Medium: case internal0Low:
        case internal1High: case internal1Medium: case internal1Low:
        case internal2High: case internal2Medium: case internal2Low:
	{
	    int dsspos = end + Constant::dss_start + 1;
	    if (end == dnalen-1) {// exon is longer than dna, right truncated exon
	      endPartProb = 1.0;  // in this case allow that there is no splice site consensus
	    } else if ((dsspos + DSS_MIDDLE - 1 < dnalen && !isPossibleDSS(dsspos)) || 
		       end + Constant::dss_start >= dnalen || 
		       orf->leftmostExonBegin(win - 1, end + Constant::dss_start, true) >= end)
		endPartProb = 0.0;
	    else 
		endPartProb = 1.0;
	    feature = seqFeatColl->getFeatureListContaining(SET_FLAG(dssF), dsspos, plusstrand);
	    if (feature) {
		while (feature) {
		    extrinsicEmiQuot *= feature->distance_faded_bonus(dsspos);
		    feature = feature->next;
		}
	    } else if (seqFeatColl->collection->hasHintsFile)
	        extrinsicEmiQuot = seqFeatColl->collection->malus(dssF) * seqFeatColl->localSSMalus(dssF, dsspos, plusstrand);
	}
	break;
        case rterminal0High: case rterminal0Medium: case rterminal0Low:
        case rterminal1High: case rterminal1Medium: case rterminal1Low:
        case rterminal2High: case rterminal2Medium: case rterminal2Low:
        case rinternal0High: case rinternal0Medium: case rinternal0Low:
        case rinternal1High: case rinternal1Medium: case rinternal1Low:
        case rinternal2High: case rinternal2Medium: case rinternal2Low:
	{
	    int asspos = end + Constant::ass_end + 1;
	    if (end == dnalen-1){ // exon is longer than dna, right truncated exon
	      endPartProb = 1.0;  // in this case allow that there is no splice site consensus
	    } else if (end + Constant::ass_end + ASS_MIDDLE < dnalen && isPossibleRASS(asspos)){
#ifdef DEBUG
		if (!seqFeatColl->validRASSPattern(sequence + asspos))
		    throw ProjectError("pattern found is not valid");
#endif
		endPartProb = 1.0;
	    } else 
		endPartProb = 0.0;
	    feature = seqFeatColl->getFeatureListContaining(SET_FLAG(assF), asspos, minusstrand);
	    if (feature)
		while (feature) {
		    extrinsicEmiQuot *= feature->distance_faded_bonus(asspos);
		    feature = feature->next;
		}
	    else if (seqFeatColl->collection->hasHintsFile)
	        extrinsicEmiQuot = seqFeatColl->collection->malus(assF) * seqFeatColl->localSSMalus(assF, asspos, minusstrand);
	}
	    break;
	default:
	    cerr << "ExonModel::viterbiAlgorithm: unknown alternative." << endl;
	    cerr << "etype: " << etype << endl;
	    endPartProb =  0.0;
    }
    return endPartProb * extrinsicEmiQuot;
}


/*
 * ===[ ExonModel::notEndPartEmiProb ]=====================================
 *
 * Probability of the emission of the exon exluding the fixed length end part.
 * This includes 
 * - the beginPart     (translation initiation, ass or reverse dss)
 * - the startSeq part (first k or less emissions after the begin part)
 * - the restSeq part  (rest of the variable length sequence emission part)
 * - the length part   (probability of the length of the biological exon)
 *
 * "right" is the last position of the inner sequence part as described in viterbiAlgorithm.
 * frameOfRight is the reading frame of the position "right".
 */

Double ExonModel::notEndPartEmiProb(int beginOfStart, int right, int frameOfRight,
				    Feature *exonparts) const {
    Double beginPartProb, restSeqProb, lenPartProb;
    Feature *feature;
    Double extrinsicQuot = 1;

    /*
     * probability of the begin part
     */
    int beginOfBioExon = beginOfStart - innerPartOffset;
    switch( etype ){
	case singleGHigh: case singleGMedium: case singleGLow:
	case initial0High: case initial0Medium: case initial0Low:
	case initial1High: case initial1Medium: case initial1Low:
	case initial2High: case initial2Medium: case initial2Low:
	    // start codon at the beginning?
	    if ((beginOfBioExon >= 0) && onStart(sequence + beginOfStart - STARTCODON_LEN)) {
		// two cases ... . the normal one with enough sequence space before the gene
		int transInitStart = beginOfBioExon - trans_init_window;
		if (transInitStart > transInitMotif->k)
		    beginPartProb = transInitMotif->seqProb(sequence+transInitStart);
		/* ... and the case where there is no place for the transInitMotif
		 * take emission probs of 1/4 for the rest up to the beginning of the seq
		 * endOfPred is negative in this case!
		 * Need this if the gene starts right after the sequence.
		 */
		else {
		    beginPartProb = pow(0.25, (double)(beginOfStart - STARTCODON_LEN));
		}
		feature = seqFeatColl->getFeatureListOvlpingRange(startF, beginOfStart-3, beginOfStart-1 , plusstrand); 
		if (feature) {
		  while (feature) {
		    if (feature->start <= beginOfStart-3 && feature->end >= beginOfStart-1)
		      extrinsicQuot *= feature->distance_faded_bonus(beginOfStart-2);
		    feature = feature->next;
		  }
		} else if (seqFeatColl->collection->hasHintsFile)
		  extrinsicQuot = seqFeatColl->collection->malus(startF);
		/*
		feature = seqFeatColl->getFeatureAt(startF, beginOfStart-1 , plusstrand);
		if (feature)
		    extrinsicQuot *= feature->bonus;
		else if (seqFeatColl->collection->hasHintsFile)
		    extrinsicQuot *= seqFeatColl->collection->malus(startF);*/
	    } else 
		beginPartProb = 0.0; 
	    break;
	case terminalHigh: case terminalMedium: case terminalLow:
	case internal0High: case internal0Medium: case internal0Low:
	case internal1High: case internal1Medium: case internal1Low:
	case internal2High: case internal2Medium: case internal2Low:
	    if (beginOfStart > 0) {
		// only a shortcut if there is no possible ass at the right position:
		if ((beginOfBioExon < 0) || 
		    ((beginOfBioExon - ASS_MIDDLE >=0 ) && !isPossibleASS(beginOfBioExon - 1)))
		    beginPartProb = 0.0;
		else 
		    beginPartProb = 1.0; // splice site is evaluated in other state
		if (beginPartProb > 0.0){
		    feature = seqFeatColl->getFeatureListContaining(SET_FLAG(assF), beginOfBioExon - 1, plusstrand);
		    if (feature) {
			while (feature) {
			    extrinsicQuot *= feature->distance_faded_bonus(beginOfBioExon - 1);
			    feature = feature->next;
			}
		    } else if (seqFeatColl->collection->hasHintsFile)
		        extrinsicQuot = seqFeatColl->collection->malus(assF) * seqFeatColl->localSSMalus(assF, beginOfBioExon - 1, plusstrand);
		}
	    } else if (beginOfStart == 0) // left truncated
	        beginPartProb = 1.0;
	    else 
	        beginPartProb = 0.0;
	    break;
	case rsingleGHigh: case rsingleGMedium: case rsingleGLow:
	case rterminal0High: case rterminal0Medium: case rterminal0Low:
	case rterminal1High: case rterminal1Medium: case rterminal1Low: 
	case rterminal2High: case rterminal2Medium: case rterminal2Low:
	    // reverse stop codon
	    if (beginOfBioExon < 0)
		beginPartProb = 0.0;
	    else if (strncmp(sequence + beginOfBioExon, "tta", 3)==0)
		beginPartProb = Constant::ochreprob;
	    else if (strncmp(sequence + beginOfBioExon, "cta", 3)==0)
		beginPartProb = Constant::amberprob;
	    else if (strncmp(sequence + beginOfBioExon, "tca", 3)==0)
		    beginPartProb = Constant::opalprob;
	    else 
		beginPartProb = 0.0;
	    // extrinsic info about a reverse stop codon?
	    if (beginPartProb > 0.0){
	       feature = seqFeatColl->getFeatureListOvlpingRange(stopF, beginOfStart-3, beginOfStart-1 , minusstrand); 
	       if (feature) {
		 while (feature) {
		   if (feature->start <= beginOfStart-3 && feature->end >= beginOfStart-1)
		     extrinsicQuot *= feature->distance_faded_bonus(beginOfStart-2);
		   feature = feature->next;
		 }
	       } else if (seqFeatColl->collection->hasHintsFile)
		 extrinsicQuot = seqFeatColl->collection->malus(stopF);
	       /*
	       feature = seqFeatColl->getFeatureAt(stopF, beginOfStart-1 , minusstrand);
		if (feature)
		    extrinsicQuot *= feature->bonus;
		else  if (seqFeatColl->collection->hasHintsFile)
		extrinsicQuot *= seqFeatColl->collection->malus(stopF);*/
	    }
	    break;
	case rinitialHigh: case rinitialMedium: case rinitialLow: 
	case rinternal0High: case rinternal0Medium: case rinternal0Low:
	case rinternal1High: case rinternal1Medium: case rinternal1Low:
	case rinternal2High: case rinternal2Medium: case rinternal2Low:
	    // reverse donor splice site
	    if (beginOfStart == 0)
	      beginPartProb = 1.0; // left truncated
	    else if (beginOfBioExon<0 || (beginOfBioExon - DSS_MIDDLE>0 && !isPossibleRDSS(beginOfBioExon - 1)))
		beginPartProb = 0.0;
	    else 
		beginPartProb = 1.0; // splice site is evaluated in other state
	    if (beginPartProb > 0.0){
		feature = seqFeatColl->getFeatureListContaining(SET_FLAG(dssF), beginOfBioExon - 1, minusstrand);
		if (feature) {
		    while (feature) {
			extrinsicQuot *= feature->distance_faded_bonus(beginOfBioExon - 1);
			feature = feature->next;
		    }
		} else if (seqFeatColl->collection->hasHintsFile)
		    extrinsicQuot = seqFeatColl->collection->malus(dssF) * seqFeatColl->localSSMalus(dssF, beginOfBioExon - 1, minusstrand);
	    }
	    break;
	default:
	    throw ExonModelError("notEndPartEmiProb: Wrong ExonType");
	    break;
    }

    if (!(beginPartProb > 0.0)) 
	return 0.0;
	
    /*
     * sequence emission probability: restSeqProb
     */

    if (beginOfStart > right) {
	/*
	 * inner sequence empty or negative length = overlapping begin and end part
	 * this is for dealing with very short exons, which would be prohibited otherwise
	 * by the length of the fixed length emission of the beginning and end part.
	 * The method here is: Use the beginning and end part even if not really applicable
	 * and correct for the fact that the emissions of the overlapping part were taken into
	 * account twice by dividing by a standard emission probability of 1/4.
	 */
	restSeqProb = POWER4TOTHE( beginOfStart - right -1);
    } else if (right-beginOfStart <= k){
	/*
	 *  inner sequence has length <= k+1, a short exon but with inner sequence of positive length
	 */
	restSeqProb = 1.0;
	Seq2Int s2i(right-beginOfStart+1);
	try {
	    if (isOnFStrand(etype)) {
	    	if (isHighState(etype)){
	    		restSeqProb = highPls[right-beginOfStart][frameOfRight][s2i(sequence + beginOfStart)];
	    	} else if(isMediumState(etype)) {
	    		restSeqProb = mediumPls[right-beginOfStart][frameOfRight][s2i(sequence + beginOfStart)];
	    	} else if (isLowState(etype)) {
	    		restSeqProb = lowPls[right-beginOfStart][frameOfRight][s2i(sequence + beginOfStart)];
	    	}
	    }
	    else
		{
	    	if (isHighState(etype)) {
	    		restSeqProb = highPls[right-beginOfStart][mod3(frameOfRight + right - beginOfStart)][s2i.rc(sequence + beginOfStart)];
	    	} else if (isMediumState(etype)) {
	    		restSeqProb = mediumPls[right-beginOfStart][mod3(frameOfRight + right - beginOfStart)]
	    			    		                                      [s2i.rc(sequence + beginOfStart)];
	    	} else if (isLowState(etype)) {
	    		restSeqProb = lowPls[right-beginOfStart][mod3(frameOfRight + right - beginOfStart)]
	    			    		                                      [s2i.rc(sequence + beginOfStart)];
	    	}

		}
	} catch (InvalidNucleotideError e) {
	    // we dont assume anything in this case, take iid uniform distribution on {a,c,g,t}
	    restSeqProb = pow(Constant::probNinCoding, right-beginOfStart + 1); // 0.25
	}
    } else {
	/*
	 * inner sequence has length > k+1, this is the normal case
	 */
	
	/*
	 * inner sequence part from endOfStart+1 to right
	 * | initial pattern | initial content model| exon content model| exon terminating model (dss)|
	 * Depending on the exon type, some of the models may not apply
	 * For short exons the precedence is inital -> terminal -> inner
	 * but terminal applies only completely or not at all.
	 */
	int endOfInitial, beginOfTerm,  endOfStart; // forward
	int endOfTerm, beginOfInitP, beginOfInitial;// reverse

	endOfStart = beginOfStart + k-1;
	beginOfInitP = right - (k-1);
	try {
	    if (k==0)
		restSeqProb = 1;
	    else
		if (isOnFStrand(etype)) // init pattern at left side
		{
			if (isHighState(etype)) {
				restSeqProb = highPls[k-1][mod3(frameOfRight-right+endOfStart)][Seq2Int(k)(sequence + beginOfStart)];
			} else if (isMediumState(etype)) {
				restSeqProb = mediumPls[k-1][mod3(frameOfRight-right+endOfStart)][Seq2Int(k)(sequence + beginOfStart)];
			} else if (isLowState(etype)) {
				restSeqProb = lowPls[k-1][mod3(frameOfRight-right+endOfStart)][Seq2Int(k)(sequence + beginOfStart)];
			}
		}
		else                  // init pattern at right side
		    {
			if (isHighState(etype)) {
				restSeqProb = highPls[k-1][mod3(frameOfRight+right-beginOfInitP)]
				[Seq2Int(k).rc(sequence + beginOfInitP)];
			} else if (isMediumState(etype)) {
				restSeqProb = mediumPls[k-1][mod3(frameOfRight+right-beginOfInitP)]
								[Seq2Int(k).rc(sequence + beginOfInitP)];
			} else if (isLowState(etype)) {
				restSeqProb = lowPls[k-1][mod3(frameOfRight+right-beginOfInitP)]
								[Seq2Int(k).rc(sequence + beginOfInitP)];
			}
		    }
	} catch (InvalidNucleotideError e) {
	    restSeqProb = pow(Constant::probNinCoding, (int) k ); // 0.25
	}

	switch( etype ){
	    case singleGHigh: case singleGMedium: case singleGLow:
		/*
		 * | initial pattern | initial content model | content model |
		 *   <---   k   --->
		 */
		endOfInitial = endOfStart + Constant::init_coding_len;
		if (endOfInitial > right)
		    endOfInitial = right;
		restSeqProb *=
		    initialSeqProb(endOfStart+1, endOfInitial, mod3(frameOfRight-right+endOfInitial)) *
		    seqProb(endOfInitial+1, right, frameOfRight);
		break;
	    case initial0High: case initial0Medium: case initial0Low: case initial1High: case initial1Medium: case initial1Low: case initial2High: case initial2Medium: case initial2Low:
		/*
		 * | initial pattern | initial content model | content model | eterminal model |
		 *   <---   k   --->   <- init_coding_len ->                  < et_coding_len >  
		 */
		endOfInitial = endOfStart + Constant::init_coding_len;
		if (endOfInitial > right) {
		    endOfInitial = right;
		    beginOfTerm = right + 1;
		} else  {
		    beginOfTerm = right - Constant::et_coding_len + 1;
		    if (beginOfTerm  <= endOfInitial )
			beginOfTerm = right + 1; // no terminal part 
		}
		restSeqProb *= 
		    initialSeqProb(endOfStart+1, endOfInitial, mod3(frameOfRight-right+endOfInitial)) *
		    seqProb(endOfInitial+1, beginOfTerm-1, mod3(frameOfRight-right+(beginOfTerm-1))) *
		    eTermSeqProb(beginOfTerm, right, frameOfRight);
		break;	
	    case internal0High: case internal0Medium: case internal0Low: case internal1High: case internal1Medium: case internal1Low: case internal2High: case internal2Medium: case internal2Low:
		/*
		 * | initial pattern | content model | eterminal model |
		 *   <---   k   --->                  < et_coding_len >  
		 */
		beginOfTerm = right - Constant::et_coding_len + 1;
		if (beginOfTerm <= endOfStart)
		    beginOfTerm = right + 1;   // exon not long enough, no terminal part
		restSeqProb *= 
		    seqProb(endOfStart+1, beginOfTerm-1, mod3(frameOfRight-right+(beginOfTerm-1))) *
		    eTermSeqProb(beginOfTerm, right, frameOfRight);
		break;
	    case terminalHigh: case terminalMedium: case terminalLow:
		/*
		 * | initial pattern | content model |
		 *   <---   k   --->
		 */
		restSeqProb *= seqProb(endOfStart+1, right, frameOfRight);
		break;
	    case rsingleGHigh: case rsingleGMedium: case rsingleGLow:
		/*
		 * | content model | initial content model | initial pattern |
		 *                   <- init_coding_len ->   <---   k   --->
		 */
		beginOfInitial = beginOfInitP - Constant::init_coding_len;
		if (beginOfInitial < beginOfStart)
		    beginOfInitial = beginOfStart;
		restSeqProb *=
		    initialSeqProb(beginOfInitial, beginOfInitP-1, mod3(frameOfRight+right-(beginOfInitP-1))) *
		    seqProb(beginOfStart, beginOfInitial-1, mod3(frameOfRight+right-(beginOfInitial-1)));
		break;
	    case rinitialHigh: case rinitialMedium: case rinitialLow:
		/*
		 * | eterminal model | content model | initial content model | initial pattern |
		 *   < et_coding_len >                 <- init_coding_len ->   <---   k   --->
		 */
		beginOfInitial = beginOfInitP - Constant::init_coding_len;
		if (beginOfInitial < beginOfStart) {
		    beginOfInitial = beginOfStart;
		    endOfTerm = beginOfStart - 1;
		} else {	    
		    endOfTerm = beginOfStart + Constant::et_coding_len - 1;
		    if (endOfTerm >= beginOfInitial)
			endOfTerm = beginOfStart - 1;   // exon not long enough, no terminal part
		}
		restSeqProb *= 
		    initialSeqProb(beginOfInitial, beginOfInitP-1, mod3(frameOfRight+right-(beginOfInitP-1))) *
		    seqProb(endOfTerm + 1, beginOfInitial-1, mod3(frameOfRight+right-(beginOfInitial-1))) *
		    eTermSeqProb(beginOfStart, endOfTerm, mod3(frameOfRight+right-endOfTerm));
		break;	
	    case rinternal0High: case rinternal0Medium: case rinternal0Low: case rinternal1High: case rinternal1Medium: case rinternal1Low: case rinternal2High: case rinternal2Medium: case rinternal2Low:
		/*
		 * | eterminal model | content model | initial pattern |
		 *   < et_coding_len >                 <---   k   ---> 
		 */
		endOfTerm = beginOfStart + Constant::et_coding_len - 1;
		if (endOfTerm >= beginOfInitP) 
		    endOfTerm = beginOfStart - 1;   // exon not long enough, no terminal part
		restSeqProb *=
		    seqProb(endOfTerm + 1, beginOfInitP-1, mod3(frameOfRight+right-(beginOfInitP-1))) *
		    eTermSeqProb(beginOfStart, endOfTerm, mod3(frameOfRight+right-endOfTerm));
		break;
	    case rterminal0High: case rterminal0Medium: case rterminal0Low: case rterminal1High: case rterminal1Medium: case rterminal1Low: case rterminal2High: case rterminal2Medium: case rterminal2Low:
		/*
		 * | content model | initial pattern | 
		 *                   <---   k   --->
		 */
		restSeqProb *=
		    seqProb(beginOfStart, beginOfInitP-1,mod3(frameOfRight+right-(beginOfInitP-1)));
		break;
	    default:
		throw ExonModelError("notEndPartEmiProb, inner sequence part: Wrong ExonType");
		break;
	}
    }
    
    /*
     * compute the probability of the length
     */
    int endOfBioExon = right + innerPartEndOffset;
    int exonLength = endOfBioExon - beginOfBioExon + 1;

    if (exonLength < 1) {
	lenPartProb = 0.0;
    } else{
	switch( etype ){
	    case singleGHigh: case singleGMedium: case singleGLow: case rsingleGHigh: case rsingleGMedium: case rsingleGLow:
		if (exonLength%3 == 0)
		{
			if (isHighState(etype)) {
				lenPartProb = 3*highlenDistSingle[exonLength];
			} else if (isMediumState(etype)) {
				lenPartProb = 3*mediumlenDistSingle[exonLength];
			} else if (isLowState(etype)) {
				lenPartProb = 3*lowlenDistSingle[exonLength];
			}
		}
		else
		    lenPartProb = 0.0;
		break;
	    case initial0High: case initial0Medium: case initial0Low: case initial1High: case initial1Medium: case initial1Low: case initial2High: case initial2Medium: case initial2Low:
		if (exonLength%3 == win && exonLength > 2)  {
			if (isHighState(etype)) {
				lenPartProb = 3*highlenDistInitial[exonLength];
			} else if (isMediumState(etype)) {
				lenPartProb = 3*mediumlenDistInitial[exonLength];
			} else if (isLowState(etype)) {
				lenPartProb = 3*lowlenDistInitial[exonLength];
			}
		}
		else
		    lenPartProb = 0.0;
		break;
	    case rinitialHigh: case rinitialMedium: case rinitialLow:
		if (exonLength > 2) {
			if (isHighState(etype)) {
				lenPartProb = 3*highlenDistInitial[exonLength];
			} else if (isMediumState(etype)) {
				lenPartProb = 3*mediumlenDistInitial[exonLength];
			} else if (isLowState(etype)) {
				lenPartProb = 3*lowlenDistInitial[exonLength];
			}
		}
		else 
		    lenPartProb = 0.0;
		break;
	    case internal0High: case internal0Medium: case internal0Low: 
	    case internal1High: case internal1Medium: case internal1Low: 
	    case internal2High: case internal2Medium: case internal2Low: 
	    case rinternal0High: case rinternal0Medium: case rinternal0Low: 
	    case rinternal1High: case rinternal1Medium: case rinternal1Low: 
	    case rinternal2High: case rinternal2Medium: case rinternal2Low:
		/*
		 * for internal and terminal exons the length and the 
		 * predecessor state must fit together, see viterbi algorithm
		 */
	    if (isHighState(etype)) {
	    	lenPartProb = 3*highlenDistInternal[exonLength];
	    } else if (isMediumState(etype)) {
	    	lenPartProb = 3*mediumlenDistInternal[exonLength];
	    } else if (isLowState(etype)) {
	    	lenPartProb = 3*lowlenDistInternal[exonLength];
	    }
		break;
	    case terminalHigh: case terminalMedium: case terminalLow:
	    {
	    	if (isHighState(etype)) {
	    		lenPartProb = 3*highlenDistTerminal[exonLength];
	    	} else if (isMediumState(etype)) {
	    		lenPartProb = 3*mediumlenDistTerminal[exonLength];
	    	} else if (isLowState(etype)) {
	    		lenPartProb = 3*lowlenDistTerminal[exonLength];
	    	}
	    }
		break;
	    case rterminal0High: case rterminal0Medium: case rterminal0Low: case rterminal1High: case rterminal1Medium: case rterminal1Low: case rterminal2High: case rterminal2Medium: case rterminal2Low:
		if (mod3(2-exonLength) == win) {
		    if (isHighState(etype)) {
		    	lenPartProb = 3*highlenDistTerminal[exonLength];
		    } else if (isMediumState(etype)) {
		    	lenPartProb = 3*mediumlenDistTerminal[exonLength];
		    } else if (isLowState(etype)) {
		    	lenPartProb = 3*lowlenDistTerminal[exonLength];
		    }
			}
		else 
		    lenPartProb = 0.0;
		break;
	    default:
		throw ExonModelError("Wrong ExonType");
		break;
	}
    }
    /*
     * Multiply a bonus/malus to extrinsicQuot for every exonpart hint that is
     * covered by this biological exon
     */
    bool exonSupport = false; // used for malus
    bool CDSSupport = false;  // used for malus
    int numEPendingInExon=0, numCPendingInExon=0, nep=0; // just used for malus. 
    bool strandOK;
    Double partBonus = 1.0;
    for (Feature *part = exonparts; part!= NULL; part = part->next){
	strandOK = part->strand == bothstrands || part->strand == STRAND_UNKNOWN || isOnFStrand(etype) == (part->strand == plusstrand);
	if (part->type == exonpartF || part->type == CDSpartF){
	    if (part->type == exonpartF && part->end >= beginOfBioExon && part->end <= endOfBioExon)
		numEPendingInExon++;
	    if (part->type == CDSpartF && part->end >= beginOfBioExon && part->end <= endOfBioExon)
		numCPendingInExon++;
	    if (strandOK){
	        if (part->start >= beginOfBioExon && part->end <= endOfBioExon){
		    partBonus *= part->bonus;
		    nep += 1;
  	        } else { // exonpart not completely contained in exon
		  if (part->type == exonpartF) {
		    if ((etype == singleGHigh || etype == singleGMedium || etype == singleGLow || 
		    	etype == initial0High || etype == initial0Medium || etype == initial0Low || 
		    	etype == initial1High || etype == initial1Medium || etype == initial1Low || 
		    	etype == initial2High || etype == initial2Medium || etype == initial2Low ||
			    etype == rsingleGHigh || etype == rsingleGMedium || etype == rsingleGLow || 
			    etype == rterminal0High || etype == rterminal0Medium || etype == rterminal0Low || 
			    etype == rterminal1High || etype == rterminal1Medium || etype == rterminal1Low || 
			    etype == rterminal2High || etype == rterminal2Medium || etype == rterminal2Low) &&
			part->end <= endOfBioExon && part->end >= beginOfBioExon){
		      partBonus *= sqrt(part->bonus);
		      nep += 1;
		    }
		    if ((etype == singleGHigh || etype == singleGMedium || etype == singleGLow || 
		    	etype == terminalHigh || etype == terminalMedium || etype == terminalLow || 
		    	etype == rsingleGHigh || etype == rsingleGMedium || etype == rsingleGLow || 
		    	etype == rinitialHigh || etype == rinitialMedium || etype == rinitialLow) &&
			part->start >= beginOfBioExon && part->start <= endOfBioExon){
		      partBonus *= sqrt(part->bonus); 
		      nep += 1;
		    }
		    // could both happen if single exon gene has overlapping exonpart at both end points
		  }
		}
	    }
	}
	/* 
	 * Multiply a bonus to extrinsicQuot for every exon/CDS hint that exactly matches these
	 * exon boundaries. Malus for unsupported exons.
	 */
	if (part->type == CDSF && part->start == beginOfBioExon && part->end == endOfBioExon && strandOK){
	    CDSSupport = true;
	    extrinsicQuot *= part->bonus;
	}
	if (part->type == exonF && strandOK) {
	    if (etype == singleGHigh || etype == singleGMedium || etype == singleGLow || etype == rsingleGHigh || etype == rsingleGMedium || etype == rsingleGLow) {
		// do nothing, since both exons ends are handled in UTR states (if there is UTR)
	    } else if (etype == internal0High || etype == internal0Medium || etype == internal0Low || 
	    	etype == internal1High || etype == internal1Medium || etype == internal1Low ||
	    	etype == internal2High || etype == internal2Medium || etype == internal2Low || 
	    	etype == rinternal0High || etype == rinternal0Medium || etype == rinternal0Low || 
	    	etype == rinternal1High || etype == rinternal1Medium || etype == rinternal1Low || 
	    	etype == rinternal2High || etype == rinternal2Medium || etype == rinternal2Low) {
		if (part->start == beginOfBioExon && part->end == endOfBioExon) {
		    exonSupport = true;
		    extrinsicQuot *= part->bonus; 
		}
	    } else if (etype == terminalHigh || etype == terminalMedium || etype == terminalLow || etype == rinitialHigh || etype == rinitialMedium || etype == rinitialLow) {
		if (part->start == beginOfBioExon && part->end > endOfBioExon) {
		    exonSupport = true;
		    extrinsicQuot *= sqrt(part->bonus); 
		}
	    } else {
		if (part->start < beginOfBioExon && part->end == endOfBioExon) {
		    exonSupport = true;
		    extrinsicQuot *= sqrt(part->bonus); 
		}
	    }
	}
    }
    extrinsicQuot *= partBonus;
    if (seqFeatColl && seqFeatColl->collection->hasHintsFile) {
	/* We have searched for extrinsic features.
	 * Then multiply the malus for each position of the exon.
	 * We should exclude those (few) positions where an exonpart ends.
	 * Exons longer than their exonpart hint have an incentive to become shorter.
	 */
      if (nep >= 5){// local malus for partially and unevenly supported CDS
	int zeroCov = seqFeatColl->numZeroCov(beginOfBioExon, endOfBioExon, CDSpartF, isOnFStrand(etype)? plusstrand : minusstrand);
	Double localPartMalus = seqFeatColl->collection->localPartMalus(CDSpartF, zeroCov, partBonus, nep);
	if (localPartMalus < 1.0/partBonus) // at least have ab initio probabilities
	  localPartMalus = 1.0/partBonus;
	//	cout << "partBonus[" << beginOfBioExon << ", " << endOfBioExon << "]= " << partBonus << " zeroCov = " << zeroCov << " localPartMalus = " << localPartMalus << endl;
	extrinsicQuot *= localPartMalus;
      }
	if (exonLength-numEPendingInExon > 0)
	    extrinsicQuot *= seqFeatColl->collection->partMalus(exonpartF, exonLength-numEPendingInExon);
	if (exonLength-numCPendingInExon > 0)
	    extrinsicQuot *= seqFeatColl->collection->partMalus(CDSpartF, exonLength-numCPendingInExon);
    }

    if (seqFeatColl && seqFeatColl->collection->hasHintsFile) {
	// We have searched for extrinsic features but have not found an exon hint.
	if (!exonSupport)
	    extrinsicQuot *= seqFeatColl->collection->malus(exonF);
	if (!CDSSupport)
	    extrinsicQuot *= seqFeatColl->collection->malus(CDSF);
    }
    return beginPartProb * restSeqProb * lenPartProb * extrinsicQuot;
}


/*
 * ===[ ExonModel::emiProbUnderModel ]=====================================
 *
 * Probability of emitting dna[begin]...dna[end] in this exon state.
 * Includes sequence emission and state length.
 */

Double ExonModel::emiProbUnderModel(int begin, int end) const {
    if (begin > end)
	return 1.0;

    /*
     * endPartProb
     */
    Double endPartProb = endPartEmiProb(end);
    
    /*
     * right is the last position of the inner sequence part
     * endOfBioExon is the last position of the biological exon state
     * beginOfStart, see viterbiAlgorithm
     */
    
    int endOfBioExon = end + baseOffset;
    int right = endOfBioExon - innerPartEndOffset;
    int beginOfStart = begin + beginPartLen;
    int beginOfBioExon = beginOfStart - innerPartOffset;
    
    if( !(endPartProb > 0) || (right < 0))
	return 0.0;
    
     /*
     * Reading frame of the position "right"
     * win=0,1,2 is the reading frame of __ endOfBioExon+1 __
     * on the reverse strand, count backwards 
     */
    int frameOfRight = isOnFStrand(etype) ? 
	mod3(win - (endOfBioExon + 1) + right) : 
	mod3(win +  endOfBioExon + 1  - right);
    
    // Determine the leftmost possible boundary of the inner sequence part.  
    int left = orf->leftmostExonBegin(frameOfRight, right, isOnFStrand(etype));
    
    if (beginOfStart < left) // e.g. a stop codon in the reading frame
	return 0.0;
    /*
     * get the extrinsic exonpart information about parts falling in this range
     * and with fitting reading frame
     */
    Feature *extrinsicexons = NULL;
    if (seqFeatColl)
	extrinsicexons = seqFeatColl->getExonListOvlpingRange(beginOfBioExon, endOfBioExon, bothstrands, mod3(right - frameOfRight));
    seqProb(-1,0,0);
    Double notEndPartProb = notEndPartEmiProb(beginOfStart, right, frameOfRight, extrinsicexons); 
    
    return endPartProb * notEndPartProb;
}


/*
 * computes the probability of the emission of the sequence from left to right
 * left and right included
 */

Double ExonModel::seqProb(int left, int right, int frameOfRight) const {
    static Double seqProb = 1;
    static int oldleft = -1, oldright = -1, oldframe = -1;
    static StateType oldtype = TYPE_UNKNOWN;

    bool reverse = !isOnFStrand(etype);
    Seq2Int s2i(k+1);
    if (left < 0) {   // new initialisation
	seqProb = 1;
        oldleft = oldright = oldframe = -1;
        oldtype = TYPE_UNKNOWN;
	return 1.0;
    }
    if (left > right) 
        return 1.0;

    if (right == oldright && frameOfRight == oldframe && left <= oldleft && etype == oldtype) {
        for (int curpos = oldleft-1; curpos >= left; curpos--){
            try {
		int f = reverse? mod3(frameOfRight+right-curpos) : mod3(frameOfRight-right+curpos);
		int pn = reverse? s2i.rc(sequence+curpos) : s2i(sequence+curpos-k);
		if (isHighState(etype)) {
			seqProb *= highemiprobs.probs[f][pn];
		} else if (isMediumState(etype)) {
			seqProb *= mediumemiprobs.probs[f][pn];
		} else if (isLowState(etype)) {
			seqProb *= lowemiprobs.probs[f][pn];
		}
            } catch (InvalidNucleotideError e) {
                seqProb *= Constant::probNinCoding; //  0.25, 1/4
            }
        }
        oldleft = left;
        return seqProb;
    }
    
    // compute everything new
    seqProb = 1.0;
    for (int curpos = right; curpos >= left; curpos--) {
	try {
	    int f = reverse? mod3(frameOfRight+right-curpos) : mod3(frameOfRight-right+curpos);
	    int pn = reverse? s2i.rc(sequence+curpos) : s2i(sequence+curpos-k);
	    if (isHighState(etype)) {
	    	seqProb  *= highemiprobs.probs[f][pn];
	    } else if(isMediumState(etype)) {
	    	seqProb  *= mediumemiprobs.probs[f][pn];
	    } else if(isLowState(etype)) {
	    	seqProb  *= lowemiprobs.probs[f][pn];
	    }
	    if (inCRFTraining && (countEnd < 0 || (curpos >= countStart && curpos <= countEnd))){
	    	if (isHighState(etype)) {
	    		highGCemiprobs[gcIdx].addCount(highGCemiprobs[gcIdx].getIndex(f,pn));
	    	} else if (isMediumState(etype)) {
	    		mediumGCemiprobs[gcIdx].addCount(mediumGCemiprobs[gcIdx].getIndex(f,pn));
	    	} else if (isLowState(etype)) {
	    		lowGCemiprobs[gcIdx].addCount(lowGCemiprobs[gcIdx].getIndex(f,pn));
	    	}
	    }
	} catch (InvalidNucleotideError e) {
	    seqProb  *= Constant::probNinCoding; // 0.25 1/4
	}
    }
    oldleft = left;
    oldright = right;
    oldframe = frameOfRight;
    oldtype = etype;
    return seqProb;
}

/*
 * eTermSeqProb
 * terminal part of exons before donor splice site
 */ 
Double ExonModel::eTermSeqProb(int left, int right, int frameOfRight) const {
//     static Double seqProb = 1; 
//     static int oldleft = -1, oldright = -1, oldframe =-1;
//     static StateType oldtype = TYPE_UNKNOWN;
    if (left > right) 
        return 1;
    bool reverse = !isOnFStrand(etype);
    if (!(right-left+1 == Constant::et_coding_len)) 
	throw ProjectError("Assertion failed in ExonModel::eTermSeqProb");
//     the first part was never called
//     if (left == oldleft && right == oldright && frameOfRight == oldframe && etype == oldtype)
// 	return seqProb;
//     else {
	//seqProb = etMotif[mod3(frameOfRight-Constant::et_coding_len + 1)]->seqProb(
	//    sequence + left);
// old method: Markov Model
	Double seqProb = 1; 
	Seq2Int s2i(k+1);
	for (int curpos = right; curpos >= left; curpos--) {
	    try {
		int f = reverse? mod3(frameOfRight+right-curpos) : mod3(frameOfRight-right+curpos);
		int pn = reverse? s2i.rc(sequence+curpos) : s2i(sequence+curpos-k);
		if (isHighState(etype)) {
			seqProb *= highetemiprobs[f][pn];
		} else if (isMediumState(etype)) {
			seqProb *= mediumetemiprobs[f][pn];
		} else if (isLowState(etype)) {
			seqProb *= lowetemiprobs[f][pn];
		}
	    } catch (InvalidNucleotideError e) {
		seqProb *= Constant::probNinCoding; // 0.25, 1/4
	    }
	}
// 	oldleft = left;
// 	oldright = right;
// 	oldframe = frame;
// 	oldtype = etype;
	return seqProb;
//     }
}

/*
 * initialSeqProb
 * first coding part of an initial or single exon (after the initial pattern)
 */ 
Double ExonModel::initialSeqProb(int left, int right, int frameOfRight) const {
    if (left > right) 
        return 1.0;
    Double seqProb = 1.0; 
    Seq2Int s2i(k+1);
    bool reverse = !isOnFStrand(etype);
    for (int curpos = right; curpos >= left; curpos--) {
	try {
	    int f = reverse? mod3(frameOfRight+right-curpos) : mod3(frameOfRight-right+curpos);
	    int pn = reverse? s2i.rc(sequence+curpos) : s2i(sequence+curpos-k);
	    if (isHighState(etype)) {
	    	seqProb *= highinitemiprobs[f][pn];
	    } else if (isMediumState(etype)) {
	    	seqProb *= mediuminitemiprobs[f][pn];
	    } else if (isLowState(etype)) {
	    	seqProb *= lowinitemiprobs[f][pn];
	    }
	} catch (InvalidNucleotideError e) {
	    seqProb *= Constant::probNinCoding; // 0.25, 1/4
	}
    }
    return seqProb;
}


/*
 * New methods
 */

bool ExonModel::isHighState(StateType type) {
	return type == singleGHigh ||
	type == initial0High ||
	type == initial1High ||
	type == initial2High ||
	type == internal0High ||
	type == internal1High ||
	type == internal2High ||
	type == terminalHigh ||
	type == rsingleGHigh ||
	type == rinitialHigh ||
	type == rinternal0High ||
	type == rinternal1High ||
	type == rinternal2High ||
	type == rterminal0High ||
	type == rterminal1High ||
	type == rterminal2High;
}

bool ExonModel::isMediumState(StateType type) {
	return type == singleGMedium ||
	type == initial0Medium ||
	type == initial1Medium ||
	type == initial2Medium ||
	type == internal0Medium ||
	type == internal1Medium ||
	type == internal2Medium ||
	type == terminalMedium ||
	type == rsingleGMedium ||
	type == rinitialMedium ||
	type == rinternal0Medium ||
	type == rinternal1Medium ||
	type == rinternal2Medium ||
	type == rterminal0Medium ||
	type == rterminal1Medium ||
	type == rterminal2Medium;
}

bool ExonModel::isLowState(StateType type) {
	return type == singleGLow ||
	type == initial0Low ||
	type == initial1Low ||
	type == initial2Low ||
	type == internal0Low ||
	type == internal1Low ||
	type == internal2Low ||
	type == terminalLow ||
	type == rsingleGLow ||
	type == rinitialLow ||
	type == rinternal0Low ||
	type == rinternal1Low ||
	type == rinternal2Low ||
	type == rterminal0Low ||
	type == rterminal1Low ||
	type == rterminal2Low;
}
