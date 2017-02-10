/**********************************************************************
 * file:    augustus.cc
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 *
 * authors: Mario Stanke, mario@gobics.de
 *
 * date    |   author      |  changes
 * --------|---------------|------------------------------------------
 * 23.01.07| Mario Stanke  | enabled redirection of cout to outfile and cerr to errfile
 **********************************************************************/

// project includes
#include "types.hh"
#include "gene.hh"
#include "genbank.hh"
#include "namgene.hh"
#include "evaluation.hh"
#include "statemodel.hh"

// standard C/C++ includes
#include <fstream>


int verbosity=1;
Boolean checkExAcc = false;
Boolean noprediction = false;
ofstream outputfile, errorfile;
string outputfilename, errorfilename;


/*
 * evaluateOnTestSet
 *
 * takes a list of annotated sequences, predicts the genes using the extrinsic information and 
 * evaluates the accuracy of the prediction
 * gene: a list of genes (annotated sequences)
 * namgene: object for the prediction algorithms
 * extrinsicFeatures: extrinsic info, is NULL if there is none
 * props: contains the parameters
 */ 
void evaluateOnTestSet(AnnoSequence *annoseq, NAMGene &namgene, FeatureCollection &extrinsicFeatures, 
		       Strand strand);

/*
 * predictOnInputSequences
 * 
 * takes a list of raw sequences and predicts the genes using the extrinsic information
 * seq: a list of sequences
 * namgene: object for the prediction algorithms
 * extrinsicFeatures: extrinsic info, is NULL if there is none
 * props: contains the parameters
 */
void predictOnInputSequences(AnnoSequence *seq, NAMGene &namgene, FeatureCollection &extrinsicFeatures, 
			     Strand strand);


/*
 * Set the global variables according to command line options or configuration file
 */

void setParameters();

/*
 * checkExtrinsicAccuracy
 *
 * takes a list of annotated sequences and determines the bonus malus of the extrinsic information
 * evaluates the bonus of the hints without cleaning for redundancies
 * the malus is not correctly computed because 
 * This is a first step to decide of which source the hints are of higher quality.
 * namgene: object for the prediction algorithms
 * extrinsicFeatures: extrinsic info
 * props: contains the parameters
 */ 
void checkExtrinsicAccuracy(AnnoSequence *annoseq, NAMGene &namgene, FeatureCollection &extrinsicFeatures);


/*
 * cutRelevantPiece
 * If predictionStart and predictionEnd are set this function cuts
 * out the piece from predictionStart to predictionEnd, and stores the offset predictionStart in the returned
 * AnnoSequence. The memory of the parameter AnnoSequence is deleted in this case, so the whole chromosome doesn't sit
 * in memory when we actually need only a small part.
 */
void cutRelevantPiece(AnnoSequence *annoseq);

/*
 * main
 */
int main( int argc, char* argv[] ){
    string     configfile;
    string     commandline;
    Strand     strand = bothstrands; // default
    int        errorcode = 0;

    LLDouble::setOutputPrecision(3);

    for (int i=0; i<argc; i++){
	commandline += argv[i];
	if (i<argc-1) 
	    commandline += " ";
    }
    
    // determination of the configuration file
    if( argc <= 1 ){
	cout << GREETING << endl << endl;
	cout << HELP_USAGE << endl;
	return 0;
    }
    try{	
	Properties::init( argc, argv );
	Constant::init();
	Gene::init();
	setParameters(); // NOTE: need Constant to be initialised first
	StateModel::init();   // set global parameters of state models

	string filename = Properties::getProperty("queryfile");
	GBProcessor gbank(filename);
	if (Gene::gff3)
	    cout << "##gff-version 3" << endl;
	    
	cout << PREAMBLE << endl;

	/*
	 * check for extrinsic information and initialise when existent
	 */
	FeatureCollection extrinsicFeatures;
   	const char *extrinsicfilename;
	try {
	  extrinsicfilename =  Properties::getProperty("hintsfile");
	} catch (...){
	  extrinsicfilename = NULL;
	  if (verbosity) 
	    cout << "# No extrinsic information on sequences given." << endl;
	}
	if (verbosity && extrinsicfilename) {
	  cout << "# reading in the file " << extrinsicfilename << " ..." << endl;
	  extrinsicFeatures.readGFFFile(extrinsicfilename);     
	  
	  if (verbosity) 
	    cout << "# Have extrinsic information about " << extrinsicFeatures.getNumSeqsWithInfo()
		 << " sequences (in the specified range). " << endl;
	}

	if (verbosity > 1) 
	    cout << "# Initialising the parameters ..." << endl;
	BaseCount::init();
	PP::initConstants();
        NAMGene namgene; // creates and initializes the states
	StateModel::readAllParameters(); // read in the parameter files: species_{igenic,exon,intron,utr}_probs.pbl

	try{
	    string strandstr = Properties::getProperty("strand");
	    if (strandstr == "forward" || strandstr == "Forward" || strandstr == "plus" || strandstr == "Plus" 
		|| strandstr == "+" || strandstr == "Watson" || strandstr == "watson" || strandstr == "w" )
		strand = plusstrand;
	    else if (strandstr == "backward" || strandstr == "Backward" || strandstr == "minus" 
		     || strandstr == "Minus" || strandstr == "-" || strandstr == "Crick" || strandstr == "crick" 
		     || strandstr == "c" || strandstr == "reverse" || strandstr == "Reverse")
		strand = minusstrand;
	    else if (strandstr == "both")
		strand = bothstrands;
	    else if (!(strandstr == ""))
		cerr << "# Unknown option for strand: " << strandstr << endl;
	} catch (...){} // take default strand
	
	if (gbank.fileType() == fasta) {
	    /*
	     * Just predict the genes for every sequence in the file.
	     */
	    if (verbosity>2) {
		if (filename == "-")
		    cout << "# Reading sequences from standard input. Assuming fasta format." << endl;
		else
		    cout << "# Looks like " << filename << " is in fasta format." << endl;
	    }
	    AnnoSequence *testsequence = gbank.getSequenceList();
	    cutRelevantPiece(testsequence);
	    predictOnInputSequences(testsequence, namgene, extrinsicFeatures, strand);
	    AnnoSequence::deleteSequence(testsequence);
	} else if (gbank.fileType() == genbank) {
	    /*
	     * Sequences were already annotated. Predict and also check the accuracy.
	     */
	    if (verbosity>2)
		cout << "# Looks like " << filename << " is in genbank format. " 
		     << "Augustus uses the annotation for evaluation of accuracy." << endl;
	    AnnoSequence *annoseq = gbank.getAnnoSequenceList();
	    cutRelevantPiece(annoseq);
	    if (!checkExAcc)
		evaluateOnTestSet(annoseq, namgene, extrinsicFeatures, strand);
	    else { // do not predict just check the accuracy of the extrinsic information
		   // without deleting for redundancies
		checkExtrinsicAccuracy(annoseq, namgene, extrinsicFeatures);
	    }
	    AnnoSequence::deleteSequence(annoseq);
	} else {
	    throw ProjectError("File format of " + filename + " not recognized.");
	}
//	if (verbosity>2)
	cout << "# command line:" << endl << "# " << commandline << endl;
    } catch( ProjectError& err ){
        cerr << "\n" <<  argv[0] << ": ERROR\n\t" << err.getMessage( ) << "\n\n";
        errorcode=1;
    } catch ( HelpException help ) {
    	cerr << help.message << endl;
    }
    if (outputfile.is_open())
	outputfile.close();
    if (errorfile.is_open())
	errorfile.close();
    return errorcode;
}

/*
 * evaluateOnTestSet
 */

void evaluateOnTestSet(AnnoSequence *annoseq, NAMGene &namgene, FeatureCollection &extrinsicFeatures, 
		       Strand strand){
    int dnaproben = 0;
    Double quotient; 
    Evaluation eval;
    Gene *genes;

    // TEMP Zeitmessung
    double total =0;

     /*// K    
     * How good was the extrinsic information.

    if (extrinsicFeatures.getNumSeqsWithInfo()){
	int numC = extrinsicFeatures.getNumCommonSeqs(annoseq);
	if (numC == 0) 
	    cout << "# WARNING: extrinsic information given but not on any of the sequences in the input set!" << endl
		 << "# Maybe different sequence names?" << endl;
	if (verbosity > 2){
	    cout << "# We have hints for " << numC << " of the sequences in the input set." << endl;
	    cout << "# Reliability of extrinsic information: " << endl;
	    extrinsicFeatures.printAccuracyForSequenceSet(annoseq, true);
	}
    }
    // to output the annotated amino acid sequence:
    if (noprediction){
	while( annoseq ){
	    //cout << ">" << annoseq->seqname << endl;
	    printGeneList(annoseq->anno->genes, annoseq, Gene::print_cds, true);
	    annoseq = annoseq->next;
	}
    } */

    while( annoseq ){
	dnaproben++;
	if (verbosity) {
	    cout << "#\n# ----- sequence number " << dnaproben << " (length = "
		 << annoseq->length << ", name = "
		 << annoseq->seqname << ") -----" << endl << "#" << endl;
	    cout << "\n# annotation: " << endl; 
	    annoseq->printGFF();
	}
	clock_t anfang, ende;
	anfang = clock();
	Evaluation eval1;

	/*
	 * check for extrinsic information about this sequence
	 */
	SequenceFeatureCollection& sfc = extrinsicFeatures.getSequenceFeatureCollection(annoseq->seqname);
	sfc.prepare(annoseq, verbosity>0 && !(Gene::gff3));
	
	bool singlestrand = false; // use not the shadow states
	try {
	     singlestrand = (Properties::getIntProperty("singlestrand") == 1);
	} catch (...) {}

	cout << "# Predicted genes for sequence number " << dnaproben <<  " on ";
	if (strand==plusstrand)
	    cout << "forward strand" << endl;
	else if (strand==minusstrand)
	    cout << "reverse strand" << endl;
	else 
	    cout << "both strands" << endl;
	if (singlestrand)
	    cout << "# Overlapping genes on opposite strand were allowed." << endl;	    
	genes = namgene.doViterbiPiecewise(sfc, annoseq, strand); 
	ende = clock();
	total += (double) (ende-anfang) / CLOCKS_PER_SEC;
	//cout << "time " << (double) (ende-anfang) / CLOCKS_PER_SEC << ", seqlen=" << annoseq->length << endl;

	eval.addToEvaluation(genes, annoseq->anno->genes, annoseq, strand, quotient);

	eval1.addToEvaluation(genes, annoseq->anno->genes, annoseq, strand, quotient);
	if (!noprediction) {
		eval1.finishEvaluation();
		eval1.printQuotients();
		eval1.print(annoseq->length);
	}
	/*
	 * clean up
	 */
	Gene::destroyGeneSequence(genes);
	if (annoseq->anno->path){
	    delete annoseq->anno->path;
	    annoseq->anno->path = NULL;           // so annoseq destructor doesn't crash
	}
	if (annoseq->anno->forwardPath){
	    delete annoseq->anno->forwardPath;
	    annoseq->anno->forwardPath = NULL;           // so annoseq destructor doesn't crash
	}	
	if (annoseq->anno->backwardPath){
	    delete annoseq->anno->backwardPath;
	    annoseq->anno->backwardPath = NULL;           // so annoseq destructor doesn't crash
	}	
	if (annoseq->anno->condensedPath); {
	    delete annoseq->anno->condensedPath;
	    annoseq->anno->condensedPath = NULL; // so annoseq destructor doesn't crash
	}
	if (annoseq->anno->condensedForwardPath); {
	    delete annoseq->anno->condensedForwardPath;
	    annoseq->anno->condensedForwardPath = NULL; // so annoseq destructor doesn't crash
	}
	if (annoseq->anno->condensedBackwardPath); {
	    delete annoseq->anno->condensedBackwardPath;
	    annoseq->anno->condensedBackwardPath = NULL; // so annoseq destructor doesn't crash
	}
	//handleViterbiVariables( namgene.getViterbiVariables() );
	annoseq = annoseq->next;
    }
    if (!noprediction) {
      eval.finishEvaluation();
      eval.printQuotients();
      eval.print();
    }
    cout << "# total time: " << total << endl;
}


/*
 * predictOnInputSequences
 */

void predictOnInputSequences(AnnoSequence *seq, NAMGene &namgene, FeatureCollection &extrinsicFeatures, 
			     Strand strand){
    int dnaproben = 0;
    int successfull = 0;
    AnnoSequence *curseq;
    //Double pathemiProb;
    Evaluation eval;
    Gene *genes;

    int numC = extrinsicFeatures.getNumCommonSeqs(seq);
    if (numC == 0 && extrinsicFeatures.getNumSeqsWithInfo() > 0) 
	cout << "# WARNING: extrinsic information given but not on any of the sequences in the input set!" << endl
	     << "# Maybe different sequence names?" << endl;
    if (verbosity>0) {
	int anz = extrinsicFeatures.getNumCommonSeqs(seq);
	cout << "# We have hints for " << anz << " sequence" << (anz!=1? "s ": " ") << "and for " 
	     << numC << " of the sequences in the input set." << endl;
    }
	    
    while( seq ){
	dnaproben++;
	curseq = seq;
	seq = seq->next;  
	curseq->next = NULL;
	if (verbosity)
	    cout << "#\n# ----- prediction on sequence number " << dnaproben << " (length = "
		 << strlen(curseq->sequence) << ", name = "
		 << curseq->seqname << ") -----" << endl << "#" << endl;
	try {
	    /*
	     * check for extrinsic information about this sequence
	     */
	    SequenceFeatureCollection& sfc = extrinsicFeatures.getSequenceFeatureCollection(curseq->seqname);
	    sfc.prepare(curseq, verbosity>0 && !(Gene::gff3));
	    bool singlestrand = false; // use not the shadow states
	    try {
		singlestrand = (Properties::getIntProperty("singlestrand") == 1);
	    } catch (...) {}
	
	    cout << "# Predicted genes for sequence number " << dnaproben <<  " on ";
	    if (strand==plusstrand)
		cout << "forward strand" << endl;
	    else if (strand==minusstrand)
		cout << "reverse strand" << endl;
	    else 
		cout << "both strands" << endl;
	    if (singlestrand)
		cout << "# Overlapping genes on opposite strand are allowed." << endl;

	    genes = namgene.doViterbiPiecewise(sfc, curseq, strand); 

	    Gene::destroyGeneSequence(genes); // don't need them anymore after they are printed
	    successfull++;
	    //handleViterbiVariables( namgene.getViterbiVariables() );
	} catch (ProjectError& err ){
	    if (successfull < 1)
		throw err;
	    else 
		cerr << "\n augustus: ERROR\n\t" << err.getMessage( ) << "\n\n";
	}
	curseq->next = seq;
    }
}


/*
 * setParameters
 */

void setParameters(){
    try { 
	outputfilename = Properties::getProperty("outfile");
    } catch (...){}
    try { 
	errorfilename = Properties::getProperty("errfile");
    } catch (...){}
    try {
	verbosity = Properties::getIntProperty("/augustus/verbosity");
    } catch (...) {} 
    try {
	checkExAcc = Properties::getBoolProperty("checkExAcc");
    } catch (ProjectError e) {
	cerr << e.getMessage();
    }
    try {
      noprediction = Properties::getBoolProperty("noprediction");
    } catch (...) {}
    
    if (outputfilename != "") {
	outputfile.open(outputfilename.c_str());
	if (outputfile){
	    streambuf * fobuf = outputfile.rdbuf();
	    // the original stdout buffer could be stored at this point and restored later
	    cout.rdbuf(fobuf);
	} else {
	    cerr << "Could not open output file " << outputfilename << ". Will use stdout instead." << endl;
	}
    }
    if (errorfilename != "") {
	errorfile.open(errorfilename.c_str());
	if (errorfile){
	    cerr.rdbuf(errorfile.rdbuf());
	} else {
	    cerr << "Could not open error file " << errorfilename << ". Will use stderr instead." << endl;
	}
    }

    if (Properties::hasProperty("translation_table"))
	GeneticCode::chooseTranslationTable(Properties::getIntProperty("translation_table"));
}


void checkExtrinsicAccuracy(AnnoSequence *annoseq, NAMGene &namgene, FeatureCollection &extrinsicFeatures){
    /*
     * How good was the extrinsic information.
     */
    if (extrinsicFeatures.getNumSeqsWithInfo()){
	if (verbosity > 1){
	    cout << "Reliability of extrinsic information: " << endl;
	    extrinsicFeatures.printAccuracyForSequenceSet(annoseq, false);
	}
    }
}


/*
 * cutRelevantPiece
 * If predictionStart and predictionEnd are set this function cuts
 * out the piece from predictionStart to predictionEnd, and stores the offset predictionStart in the returned
 * AnnoSequence. The memory of the complete sequence of AnnoSequence is deleted in this case, so the whole chromosome doesn't sit
 * in memory when we actually need only a small part.
 */
void cutRelevantPiece(AnnoSequence *annoseq){
    int predictionStart, predictionEnd;
    int seqlen = annoseq->length;
    try {
	predictionStart = Properties::getIntProperty( "predictionStart" ) - 1;
    } catch (...) {
	predictionStart = 0;
    }
    try {
	predictionEnd = Properties::getIntProperty( "predictionEnd" ) - 1;
    } catch (...) {
	predictionEnd = seqlen - 1;
    }

    if (predictionStart != 0 || predictionEnd != seqlen - 1) {
	if (predictionStart < 0)
	    predictionStart = 0;
	if (predictionEnd > seqlen - 1)
	    predictionEnd = seqlen -1;
	if (predictionStart >= seqlen) 
	    throw ProjectError("predictionStart (" + itoa(predictionStart + 1) 
			       + ") is larger than sequence length (" + itoa(seqlen) 
			       + "). No predictions made.");
	if (predictionEnd < predictionStart)
	    throw ProjectError("predictionEnd (" + itoa(predictionEnd + 1) + 
			       ") is smaller than predictionStart (" + itoa(predictionStart + 1)
			       + "). No predictions made!");
	if (annoseq->next) {
	    cerr << "Warning: predictionStart or predictionEnd set but input consists of more than one sequence." << endl
		 << "Prediction will be made only on first sequence." << endl;
	    annoseq->next = 0;
	}

	annoseq->length = predictionEnd - predictionStart + 1;
	char *seq = newstrcpy(annoseq->sequence + predictionStart, annoseq->length);
	delete [] annoseq->sequence;
	annoseq->sequence = seq;
	annoseq->offset = predictionStart;
    }
}
