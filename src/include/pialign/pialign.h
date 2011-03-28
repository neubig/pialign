#ifndef PIALIGN_H__
#define PIALIGN_H__

#include "pydist.h"
#include "dirichletdist.h"
#include "base-measure.h"
#include "model-base.h"
#include "parse-chart.h"
#include "gng/string.h"
#include "gng/symbol-set.h"
#include <tr1/unordered_map>
#include <algorithm>
#include <fstream>
#include <iostream>

namespace pialign {

class PIAlign {

protected:

    // adjustable parameters
    int samples_;           // the number of samples to take
    int sampRate_;          // the distance between samples  
    int burnIn_;            // the number of samples to throw out 
    
    int babySteps_;         // the number of baby steps to take 
    int babyStepLen_;    // the length of each baby step 
    int annealSteps_;       // the number of baby steps to take 
    int annealStepLen_;  // the length of each baby step 
    int wordIters_;       // the number number of iterations to do with the word-based model (0)
    Prob annealLevel_;      // the current level of annealing
    int maxPhraseLen_;     // the maximum phrase size  
    int printMax_;         // the maximum phrase size to print  
    int printMin_;         // the minimum phrase size to print  
    int maxSentLen_;       // the maximum size of a sentence  
    int histWidth_;        // the width of the histogram to use 
    Prob probWidth_;       // the width of the probability beam to use 
    int modelType_;        // which model to use
    static const int MODEL_HIER = 0;
    static const int MODEL_FLAT = 1;
    static const int MODEL_LENGTH = 2;
    bool forceWord_;        // force the aligner to generate at-most-one alignments (phrase alignments according to the model are bracketed with {}) 
    Prob avgPhraseLen_;        // the average length of an E phrase (for flat models, default 1)
    Prob nullProb_;            // the probability null values
    Prob termStrength_, termPrior_; // the strength of the type dist and prior of the terminal
    Prob defDisc_, defStren_; // the default discounts and strengths of the Pitman-Yor process
    
    int baseType_;             // which base measure to use (default BASE_UNI)
    static const int BASE_UNI = 0;
    static const int BASE_MODEL1 = 1;
    static const int BASE_MODEL1G = 2;
    bool monotonic_;          // disable inversions

    // input/output parameters
    const char* eFile_;     // the english (target language) file
    const char* fFile_;     // the foreign (source language) file
    const char* prefix_;    // the prefix to use for output
    const char* le2fFile_;      // the file for e2f probabilities
    const char* lf2eFile_;     // the file for f2e probabilities

    // buffers and constants
    bool onSample_; // whether we are currently calculating a sample or not
    std::vector<Prob> patternBuffer_; // a buffer holding the log probabilities of each pattern type
    std::vector<Prob> poisProbs_; // a buffer holding the log probabilities of the Poisson dist
    ParseChart chartTemp_; // the chart of the current probabilities

                           //  the base distribution
    std::ostream *sampleOut_, *phraseOut_; // the file to write to
    // int currELen_, currFLen_; // the current lengths of e and f
    int currEMultiplier_; // the current lengths of e and f

    // corpora 
    // English and French corpora
    Corpus eCorpus_, fCorpus_;
    // corpus of the joint pairs used in each sentence
    Corpus aCorpus_;
    // corpora of the types of non-terminals and heads
    std::vector<int> tCorpus_, hCorpus_;

    // vocabs
    WordSymbolSet eVocab_, fVocab_;
    StringWordMap ePhrases_, fPhrases_;
    PairWordMap jointPhrases_;

    // models
    BaseMeasure* base_;
    ProbModel* model_;

    // average derivations for each phrase pair
    //  (used in printing the phrase table)
    std::vector<Prob> derivations_;

    // information variables    
    Prob likelihood_; // the likelihood of the current pass
    // time taken on each step
    double timeRemove_, timeInit_, timeBase_, timeGen_, timeFor_, timeSamp_, timeAll_;
    int totalBeam_, totalBeamTimes_;

    void trimPhraseDic(StringWordMap & dic, std::vector<WordId> & idMap);
    void trim();

public:


    PIAlign() : samples_(1), sampRate_(1), burnIn_(9),
        babySteps_(1), babyStepLen_(0), annealSteps_(1), annealStepLen_(0), wordIters_(0),
        maxPhraseLen_(3), printMax_(7), printMin_(1), maxSentLen_(40), histWidth_(0), probWidth_(1e-10),
        modelType_(MODEL_HIER), forceWord_(true), avgPhraseLen_(0.01), nullProb_(0.01), 
        termStrength_(1), termPrior_(1.0/3.0),
        defDisc_(-1), defStren_(-1),
        baseType_(BASE_MODEL1G), monotonic_(false),
        eFile_(0), fFile_(0), prefix_(0), le2fFile_(0), lf2eFile_(0), 
        patternBuffer_(3), sampleOut_(0),
        base_(0), model_(0), derivations_()
        { }

    ~PIAlign() {
        delete base_; delete model_;
    }

    ////////////////////
    // initialization //
    ////////////////////

    // load a corpus
    Corpus loadCorpus(std::string file, gng::SymbolSet< std::string, WordId > & vocab, WordId boost);

    // load the configuration
    void loadConfig(int argc, const char** argv);

    // load the corpora
    void loadCorpora();

    // allocate memory and train word-based probabilities
    void initialize();
    void loadModelOne(const char* file, WordSymbolSet & eVocab, WordSymbolSet & fVocab, PairProbMap & modelT,bool forward);
    void trainModelOne(const Corpus & es, const Corpus & fs, int eSize, int fSize, PairProbMap & modelT);

    // do phrase training
    void train();

    // bases
    void addHierBases(const WordString & e, const WordString & f);
    void addFlatNulls(const WordString & e, const WordString & f, bool forward);
    std::vector<Prob> spanModelOne(const WordString & e, const WordString & f);
    void addFlatBases(const WordString & e, const WordString & f);

    // find all edges in a sentence
    std::vector<LabeledEdge> findEdges(const WordString & str, const StringWordMap & dict);

    // generative
    void addGenerativeProbs(const WordString & e, const WordString & f, ParseChart & chart, SpanProbMap & genChart);

    // add forward probabilities
    void addForwardProbs(int eLen, int fLen, ParseChart & chart);

    // do the actual sampling
    SpanNode * sampleTree(int sent, const Span & mySpan, WordString & sentIds, int* typeIds, const ParseChart & chart, const SpanProbMap & genChart, const SpanProbMap & baseChart, bool add);
    // WordId sampleTree(int sent, const Span & mySpan, WordString & ids, int* tCounts, const ParseChart & chart, const SpanProbMap & genChart, const SpanProbMap & baseChart, bool add);
    
    // print a span
    void printSpan(const WordString & e, const WordString & f, const Span & mySpan, std::ostream & out, const char* phraseSep = " ||| ", const char* wordSep = " ", const char* phraseBeg = "((( ", const char* phraseEnd = " )))");
    std::string printSpan(const WordString & e, const WordString & f, const Span & mySpan, const char* phraseSep = " ||| ", const char* wordSep = " ",  const char* phraseBeg = "((( ", const char* phraseEnd = " )))");

    // *** sample algorithms
    void removeSample(int s);
    SpanNode * buildSample(int s, ParseChart & chart);
    WordId addSample(const WordString & e, const WordString & f, const SpanNode * myNode, int* tCounts, WordString & sentIds);
    void printSample(const WordString & e, const WordString & f, const SpanNode * myNode);

    // print the phrase table
    void printPhraseTable(std::ostream & os);

    /////////////
    // getters //
    /////////////

    // vocab
    WordSymbolSet & getEVocab() { return eVocab_; }
    const WordSymbolSet & getEVocab() const { return eVocab_; }
    void setEVocab(const WordSymbolSet & v) { eVocab_ = v; }
    WordSymbolSet & getFVocab() { return fVocab_; }
    const WordSymbolSet & getFVocab() const { return fVocab_; }
    void setFVocab(const WordSymbolSet & v) { fVocab_ = v; }

    // corpus
    const WordString & getESentence(int i) const { return eCorpus_[i]; }
    const Corpus & getECorpus() const { return eCorpus_; }
    void setECorpus(const Corpus & v) { eCorpus_ = v; }
    int getCorpusSize() const { return eCorpus_.size(); }
    const WordString & getFSentence(int i) const { return fCorpus_[i]; }
    const Corpus & getFCorpus() const { return fCorpus_; }
    void setFCorpus(const Corpus & v) { fCorpus_ = v; }

};

}

#endif
