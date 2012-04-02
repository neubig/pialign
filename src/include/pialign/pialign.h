#ifndef PIALIGN_H__
#define PIALIGN_H__

#include "pialign/pydist.h"
#include "pialign/dirichletdist.h"
#include "pialign/base-measure.h"
#include "pialign/model-base.h"
#include "pialign/parse-chart.h"
#include "pialign/look-base.h"
#include "gng/string.h"
#include "gng/symbol-set.h"
#include <tr1/unordered_map>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <pthread.h>

namespace pialign {

class PIAlign;

// information about the current iteration
class JobDetails {
public:
    Prob likelihood;
    int words,sentences,accepted;
    int totalBeam, totalBeamTimes;
    std::vector<int> beamWidths; // DEBUG
    double timeInit, timeBase, timeGen, timeLook, timeFor, timeSamp, timeRemove, timeAll;
    Prob oldProp, newProp, chartProb;
    void reset() {
        likelihood = 0; words = 0; sentences=0; accepted=0;
        totalBeam = 0; totalBeamTimes = 0;
        timeRemove = 0; timeInit = 0; timeBase = 0;
        timeGen = 0; timeLook = 0; timeFor = 0; timeSamp = 0; timeAll = 0;
        oldProp = 0; newProp = 0; chartProb = 0;
        fill(beamWidths.begin(), beamWidths.end(), 0);
    }
    JobDetails & operator+=(const JobDetails &rhs) {
        likelihood += rhs.likelihood; words += rhs.words;
        sentences += rhs.sentences; accepted += rhs.accepted;
        totalBeam += rhs.totalBeam; totalBeamTimes += rhs.totalBeamTimes;
        timeRemove += rhs.timeRemove; timeInit += rhs.timeInit; 
        timeBase += rhs.timeBase; timeGen += rhs.timeGen;
        timeLook += rhs.timeLook; 
        timeFor += rhs.timeFor; timeSamp += rhs.timeSamp; 
        timeAll += rhs.timeAll;
        oldProp += rhs.oldProp; newProp += rhs.newProp;
        chartProb += rhs.chartProb;
        if(beamWidths.size() < rhs.beamWidths.size())
            beamWidths.resize(rhs.beamWidths.size());
        for(unsigned i = 0; i < rhs.beamWidths.size(); i++)
            beamWidths[i] += rhs.beamWidths[i];
        return *this;
    }
    void printStats(std::ostream & out) {
        out << " Likelihood="<< likelihood/words <<std::endl;       
        out << " Time="<<timeAll<<"s (r="<<timeRemove<<", i="<<timeInit<<", b="<<timeBase<<", g="<<timeGen<<", l="<<timeLook<<", f="<<timeFor<<", s="<<timeSamp<<")"<<std::endl;
        out << " Avg. Beam="<<(double)totalBeam/totalBeamTimes<<", Accepted="<<((double)accepted/sentences*100)<<"%"<<std::endl;
        // out << " Beam widths:";
        // for(unsigned i = 1; i < beamWidths.size(); i++)
        //     out << " "<<i<<"="<<(double)beamWidths[i]/sentences;
        // out << std::endl;
    }
};

// a data structure representing a single job
class BuildJob {
public:
    ParseChart chart;
    LookAhead* lookAhead;
    std::vector<int>::iterator begin, end;
    std::vector<SpanNode*>::iterator beginOld, beginNew;
    pthread_t thread;
    PIAlign* pialign;
    JobDetails details;
    
    BuildJob() : lookAhead(0) { }
    ~BuildJob() {
        if(lookAhead)
            delete lookAhead;
    }

};

class PIAlign {

protected:

    // adjustable parameters
    int samples_;           // the number of samples to take
    int sampRate_;          // the distance between samples  
    int burnIn_;            // the number of samples to throw out 
    
    int babySteps_;        // the number of baby steps to take 
    int babyStepLen_;      // the length of each baby step 
    int annealSteps_;      // the number of baby steps to take 
    int annealStepLen_;    // the length of each baby step 
    int batchLen_;         // the length of each batch step
    int numThreads_;       // the number of threads to use
    bool shuffle_;         // shuffle the order of sentences for sampling
    int wordIters_;        // the number number of iterations to do with the word-based model (0)
    Prob annealLevel_;     // the current level of annealing
    int maxPhraseLen_;     // the maximum phrase size  
    int printMax_;         // the maximum phrase size to print  
    int printMin_;         // the minimum phrase size to print  
    int maxSentLen_;       // the maximum size of a sentence  
    Prob probWidth_;       // the width of the probability beam to use 
    bool useQueue_;        // whether to use a queue or exhaustive search
    int lookType_;             // which look-ahead to use (default LOOK_NONE)
    static const int LOOK_NONE = 0;
    static const int LOOK_IND = 1;
    static const int LOOK_INDADD = 2;

    int modelType_;        // which model to use
    static const int MODEL_HIER = 0;
    static const int MODEL_FLAT = 1;
    static const int MODEL_LENGTH = 2;
    bool forceWord_;        // force the aligner to generate at-most-one alignments (phrase alignments according to the model are bracketed with {}) 
    Prob avgPhraseLen_;        // the average length of an E phrase (for flat models, default 1)
    Prob nullProb_;            // the probability null values
    Prob termStrength_, termPrior_; // the strength of the type dist and prior of the terminal
    Prob defDisc_, defStren_; // the default discounts and strengths of the Pitman-Yor process
    
    int baseType_;             // which base measure to use (default BASE_MODEL1G)
    static const int BASE_UNI = 0;
    static const int BASE_MODEL1 = 1;
    static const int BASE_MODEL1G = 2;
    static const int BASE_PHRASECOOC_LL = 3;
    bool monotonic_;          // disable inversions
    bool doReject_;           // whether or not to do MH
    Prob coocDisc_;           // the amount to discount the coocurrence

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
    std::vector<BuildJob> buildJobs_; // sample building jobs

    // the base distribution
    int currEMultiplier_; // the current lengths of e and f

    // corpora 
    // English and French corpora
    Corpus eCorpus_, fCorpus_;
    // a corpus of derivation trees
    std::vector<SpanNode*> nCorpus_;

    // vocabs
    WordSymbolSet eVocab_, fVocab_;
    StringWordSet ePhrases_, fPhrases_;
    PairWordSet jointPhrases_;

    // models
    BaseMeasure* base_;
    ProbModel* model_;

    // average derivations for each phrase pair
    //  (used in printing the phrase table)
    std::vector<Prob> derivations_;

    void trimPhraseDic(StringWordSet & dic, std::vector<WordId> & idMap);
    void trim();

public:


    PIAlign() : samples_(1), sampRate_(1), burnIn_(9),
        babySteps_(1), babyStepLen_(0), annealSteps_(1), annealStepLen_(0), 
        batchLen_(1), numThreads_(1), shuffle_(true), wordIters_(0),
        maxPhraseLen_(3), printMax_(7), printMin_(1), maxSentLen_(40), probWidth_(1e-4), useQueue_(false), lookType_(LOOK_IND),
        modelType_(MODEL_HIER), forceWord_(true), avgPhraseLen_(0.01), nullProb_(0.01), 
        termStrength_(1), termPrior_(1.0/3.0),
        defDisc_(-1), defStren_(-1),
        baseType_(BASE_MODEL1G), monotonic_(false), doReject_(false),
        coocDisc_(1.0), 
        eFile_(0), fFile_(0), prefix_(0), le2fFile_(0), 
        lf2eFile_(0), patternBuffer_(3),
        base_(0), model_(0), derivations_()
        { }

    ~PIAlign() {
        delete base_; delete model_;
    }

    ////////////////////
    // initialization //
    ////////////////////

    // load a corpus
    void loadCorpus(Corpus & ret, std::string file, gng::SymbolSet< std::string, WordId > & vocab, WordId boost);

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

    // generative
    void addGenerativeProbs(const WordString & e, const WordString & f, ParseChart & chart, SpanProbMap & genChart) const;

    // add forward probabilities
    void addForwardProbs(int eLen, int fLen, ParseChart & chart, const SpanSet & preserve, const LookAhead & look, Prob pWidth, JobDetails & jd) const;

    // do the actual sampling
    std::pair<SpanNode*,Prob> sampleTree(int sent, const Span & mySpan, const ParseChart & chart, const SpanProbMap & genChart, const SpanProbMap & baseChart, bool add, SpanNode* actNode) const;
    
    // print a span
    void printSpan(const WordString & e, const WordString & f, const Span & mySpan, std::ostream & out, const char* phraseSep = " ||| ", const char* wordSep = " ", const char* phraseBeg = "((( ", const char* phraseEnd = " )))") const;
    std::string printSpan(const WordString & e, const WordString & f, const Span & mySpan, const char* phraseSep = " ||| ", const char* wordSep = " ",  const char* phraseBeg = "((( ", const char* phraseEnd = " )))") const;

    // *** sample algorithms
    // void *buildSamples(void* ptr);
    SpanNode * buildSample(int s, ParseChart & chart, LookAhead * lookAhead, Prob pWidth, JobDetails & jd, SpanNode* actNode) const;
    // WordId addSample(const WordString & e, const WordString & f, SpanNode * myNode);
    void printSample(const WordString & e, const WordString & f, const SpanNode * myNode, std::ostream & sampleOut, bool debug);

    void buildSpans(SpanNode* node);
    void moveRight(SpanNode* node, int e, int f);

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
    Prob getProbWidth() const { return probWidth_; }

};

}

#endif
