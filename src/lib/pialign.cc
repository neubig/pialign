

// #include "statecollection.h"
#include <fstream>
#include <cmath>
#include <sys/time.h>

#include "pialign/pialign.h"
#include "pialign/base-model1.h"
#include "pialign/base-unigram.h"
#include "pialign/base-phrasecooc.h"

#include "pialign/model-hier.h"
#include "pialign/model-flat.h"
#include "pialign/model-length.h"

#include "pialign/look-none.h"
#include "pialign/look-ind.h"

#ifdef COMPRESS_ON
#include "pialign/compress_stream.hpp"
#endif

#ifdef HAVE_CONFIG_H 
#include "pialign/config.h"
#endif

using namespace std;
using namespace std::tr1;
using namespace pialign;
using namespace gng;

void dieOnHelp(const char* err) {
#ifdef HAVE_CONFIG_H
    cerr << "---pialign ver. " << PACKAGE_VERSION "---" << endl;
#else
    cerr << "---pialign---" << endl;
#endif
cerr << " A tool for unsupervised Bayesian alignment using phrase-based ITGs" << endl
<< "  By Graham Neubig" << endl << endl
<< "Options:" << endl << endl
<< "~~~ Input/Output ~~~" << endl
<< "" << endl
<< "Usage: pialign [OPTIONS] FFILE EFILE PREFIX" << endl
<< "" << endl
<< " FFILE is the foreign input corpus" << endl
<< " EFILE is the english input corpus" << endl
<< " PREFIX is the prefix that will be used for the output" << endl
<< "" << endl
<< "Other input:" << endl
<< " -le2f         A file containing the lexicon probabilities for e2f" << endl
<< " -lf2e         A file containing the lexicon probabilities for f2e" << endl
<< "               (These can be used with \"-base m1\" or \"-base m1g\" but are not necessary)" << endl
<< "" << endl
<< "~~~ Model Parameters ~~~" << endl
<< "" << endl
<< " -model        Model type (hier/len/flat, default: hier)" << endl
<< "" << endl
<< " -avgphraselen A parameter indicating the expected length of a phrase." << endl
<< "               default is small (0.01) to prevent overly long alignments" << endl
<< " -base         The type of base measure to use (m1g is generally best)." << endl
<< "               'm1g'=geometric mean of model 1, 'm1'=arithmetic mean of model 1," << endl
<< "               'uni'=simple unigrams, 'coocll'=log-linear interpolation of phrase" << endl
<< "               cooccurrence probabilities in both directions (default 'm1g')" << endl
<< " -coocdisc     How much to discount the cooccurrence for phrasal base measures" << endl
<< " -cooccut      Cut off coocurrence probabilities at this level" << endl
<< " -defstren     Fixed strength of the PY process (default none)" << endl
<< " -defdisc      Fixed discount of the PY process (default none)" << endl
<< " -nullprob     The probability of a null alignment (default 0.01)" << endl
<< " -noremnull    Do not remember nulls in the phrase table" << endl
<< " -termprior    The prior probability of generating a terminal (0.33)" << endl
<< " -termstren    Strength of the type distribution (default 1)" << endl
#ifdef MONOTONIC_ON
<< " -monotonic    Do not allow reordering" << endl
#endif
<< "" << endl
<< "~~~ Phrase Table ~~~" << endl
<< "" << endl
<< " -maxphraselen The maximum length of a minimal phrase (default 7)" << endl
<< " -maxsentlen   The maximum length of sentences to use (default 40)" << endl
<< " -printmax     The maximum length of phrases included in the phrase table (default 7)" << endl
<< " -printmin     The minimal length of phrases included in the phrase table (default 1)" << endl
<< " -noword       Output only phrase alignments (do not force output of word alignments)" << endl
<< "" << endl
<< "~~~ Inference Parameters ~~~" << endl
<< "" << endl
<< " -burnin       The number of burn-in iterations (default 9)" << endl
<< " -histwidth    The width of the histogram pruning to use (default none)" << endl
<< " -probwidth    The width of the probability beam to use (default 1e-10)" << endl
<< " -lookahead    The type of lookahead function to use:" << endl
<< "               'none'=no look-ahead, 'ind'=independently calculate both sides" << endl
<< " -samps        The number of samples to take (default 1)" << endl
<< " -samprate     Take samples every samprate turns (default 1)" << endl
<< " -worditers    The number of iterations to perform with a word-based model (default 0)" << endl
<< " -noshuffle    Don't shuffle the order of the sentences" << endl
<< " -batchlen     The number of sentences to process in a single batch" << endl
<< " -threads      The number of threads to use (must be <= -batchlen)" << endl << endl;
    if(err)
        cerr << endl << "Error: " << err << endl;
    exit(1);
}

void PIAlign::loadConfig(int argc, const char** argv) {
    int i;
    bool rememberNull = true;
    for(i = 1; i < argc; i++) {
        if(argv[i][0] == '-') {
            if(!strcmp(argv[i],"-samps"))               samples_ = atoi(argv[++i]);
            else if(!strcmp(argv[i],"-maxphraselen"))   maxPhraseLen_ = atoi(argv[++i]);
            else if(!strcmp(argv[i],"-printmax"))       printMax_ = atoi(argv[++i]);
            else if(!strcmp(argv[i],"-printmin"))       printMin_ = atoi(argv[++i]);
            else if(!strcmp(argv[i],"-avgphraselen"))   avgPhraseLen_ = atof(argv[++i]);
            else if(!strcmp(argv[i],"-maxsentlen"))     maxSentLen_ = atoi(argv[++i]);
            else if(!strcmp(argv[i],"-burnin"))         burnIn_ = atoi(argv[++i]);
            else if(!strcmp(argv[i],"-samprate"))       sampRate_ = atoi(argv[++i]);
            else if(!strcmp(argv[i],"-histwidth"))      histWidth_ = atoi(argv[++i]);
            else if(!strcmp(argv[i],"-probwidth"))      probWidth_ = atof(argv[++i]);
            else if(!strcmp(argv[i],"-defdisc"))        defDisc_ = atof(argv[++i]);
            else if(!strcmp(argv[i],"-defstren"))       defStren_ = atof(argv[++i]);
            else if(!strcmp(argv[i],"-noshuffle"))      shuffle_ = false;
            else if(!strcmp(argv[i],"-noword"))         forceWord_ = false;
            else if(!strcmp(argv[i],"-noremnull"))      rememberNull = false;
            else if(!strcmp(argv[i],"-monotonic"))      monotonic_ = true;
            else if(!strcmp(argv[i],"-babysteps"))      babySteps_ = atoi(argv[++i]);
            else if(!strcmp(argv[i],"-babysteplen"))    babyStepLen_ = atoi(argv[++i]);
            else if(!strcmp(argv[i],"-annealsteps"))    annealSteps_ = atoi(argv[++i]);
            else if(!strcmp(argv[i],"-annealsteplen"))  annealStepLen_ = atoi(argv[++i]);
            else if(!strcmp(argv[i],"-batchlen"))       batchLen_ = atoi(argv[++i]);
            else if(!strcmp(argv[i],"-threads"))        numThreads_ = atoi(argv[++i]);
            else if(!strcmp(argv[i],"-worditers"))      wordIters_ = atoi(argv[++i]);
            else if(!strcmp(argv[i],"-le2f"))           le2fFile_ = argv[++i];
            else if(!strcmp(argv[i],"-lf2e"))           lf2eFile_ = argv[++i];
            else if(!strcmp(argv[i],"-model")) {
                ++i;
                if(!strcmp(argv[i],"hier")) modelType_ = MODEL_HIER;
                else if(!strcmp(argv[i],"flat")) modelType_ = MODEL_FLAT;
                else if(!strcmp(argv[i],"len")) modelType_ = MODEL_LENGTH;
                else {
                    ostringstream oss;
                    oss << "Unknown model argument "<<argv[i];
                    dieOnHelp(oss.str().c_str());
                }
            }
            else if(!strcmp(argv[i],"-base")) {
                ++i;
                if(!strcmp(argv[i],"m1")) baseType_ = BASE_MODEL1;
                else if(!strcmp(argv[i],"m1g")) baseType_ = BASE_MODEL1G;
                else if(!strcmp(argv[i],"uni")) baseType_ = BASE_UNI;
                else if(!strcmp(argv[i],"coocll")) baseType_ = BASE_PHRASECOOC_LL;
                else {
                    ostringstream oss;
                    oss << "Unknown base argument "<<argv[i];
                    dieOnHelp(oss.str().c_str());
                }
            }
            else if(!strcmp(argv[i],"-lookahead")) {
                ++i;
                if(!strcmp(argv[i],"none")) lookType_ = LOOK_NONE;
                else if(!strcmp(argv[i],"ind")) lookType_ = LOOK_IND;
                else if(!strcmp(argv[i],"indadd")) lookType_ = LOOK_INDADD;
                else {
                    ostringstream oss;
                    oss << "Unknown base argument "<<argv[i];
                    dieOnHelp(oss.str().c_str());
                }
            }
            else if(!strcmp(argv[i],"-coocdisc"))       coocDisc_ = atof(argv[++i]);
            else if(!strcmp(argv[i],"-cooccut"))        coocCut_ = atof(argv[++i]);
            else if(!strcmp(argv[i],"-nullprob"))       nullProb_ = atof(argv[++i]);
            else if(!strcmp(argv[i],"-termstren"))      termStrength_ = atof(argv[++i]);
            else if(!strcmp(argv[i],"-termprior"))      termPrior_ = atof(argv[++i]);
            else {
                ostringstream oss;
                oss << "Unknown argument "<<argv[i];
                dieOnHelp(oss.str().c_str());
            }
        }
        else if(!eFile_)
            eFile_ = argv[i];
        else if(!fFile_)
            fFile_ = argv[i];
        else if(!prefix_)
            prefix_ = argv[i];
        else
            dieOnHelp("Too many arguments!");
    }

    if(modelType_ == MODEL_HIER) model_ = new HierModel();
    else if(modelType_ == MODEL_FLAT) model_ = new FlatModel();
    else if(modelType_ == MODEL_LENGTH) model_ = new LengthModel();
    model_->setRememberNull(rememberNull);
    model_->setTerminalStrength(termStrength_);
    model_->setTerminalPrior(termPrior_);

    // sanity check
    if(!prefix_)
        dieOnHelp("Must specify the foreign file, english file, and output prefix");
    if(probWidth_ < 0 || probWidth_ >= 1)
        dieOnHelp("Probability width must be between zero and one!");
    if(nullProb_ <= 0 || nullProb_ >= 1)
        dieOnHelp("Null probability must be between zero and one!");
#ifndef MONOTONIC_ON
    if(monotonic_)
        dieOnHelp("pialign was not compiled with support for monotone decoding");
#endif
    probWidth_ = log(probWidth_);

    if(batchLen_ < 1)
        dieOnHelp("Batch length must be at least 1");
    if(numThreads_ > batchLen_)
        dieOnHelp("The number of threads must be <= the batch length");

}

double timeDifference(const timeval & s, const timeval & e) {
    return (e.tv_sec-s.tv_sec)+(e.tv_usec-s.tv_usec)/1000000.0;
}

// load the corpus
void PIAlign::loadCorpus(Corpus & ret, string file, SymbolSet< string, WordId > & vocab, WordId boost) {
    vocab.getId("",true);
    ifstream ifs(file.c_str());
    if(!ifs) {
        ostringstream oss;
        oss << "IO Error: Couldn't open "<<file;
        throw runtime_error(oss.str());
    }
    string line,buff;
    while(getline(ifs, line)) {
        istringstream iss(line);
        ret.startSentence();
        while(iss >> buff)
            ret.addWord(vocab.getId(buff,true)+boost);
        ret.endSentence();
    }
    ret.makeSentences();
}

void PIAlign::loadCorpora() {

    // load the corpora
    loadCorpus(eCorpus_, eFile_, eVocab_, 0);
    loadCorpus(fCorpus_, fFile_, fVocab_, eVocab_.size());
    // remove sentences that are too long
    int total = 0, maxLen = 0;
    for(unsigned i = 0; i < eCorpus_.size(); i++) {
        int myLen = max(eCorpus_[i].length(),fCorpus_[i].length());
        if(myLen > maxSentLen_) {
            eCorpus_[i].setLength(0);
            fCorpus_[i].setLength(0);
        } else {
            total++;
            maxLen = max(myLen,maxLen);
        }
    }
    if(eCorpus_.size() != fCorpus_.size())
       throw std::runtime_error("Corpus sizes don't match");
    cerr << "Loaded corpus: "<<total<<"/"<<eCorpus_.size()<<" sentences used"<<endl;
    
    // allocate other corpora
    nCorpus_ = vector<SpanNode*>(eCorpus_.size(),0);
    
    // initialize the charts
    buildJobs_.resize(numThreads_);
    for(int i = 0; i < numThreads_; i++) {
        buildJobs_[i].chart.initialize(maxLen,maxLen);
        buildJobs_[i].pialign = this;
        // -------------------- make the look-ahead -----------------------
        // make the look-ahead
        if(lookType_ == LOOK_NONE)
            buildJobs_[i].lookAhead = new LookAheadNone();
        else if(lookType_ == LOOK_IND || lookType_ == LOOK_INDADD) {
            buildJobs_[i].lookAhead = new LookAheadInd();
            if(lookType_ == LOOK_INDADD)
                ((LookAheadInd*)buildJobs_[i].lookAhead)->setAdd(true);
        }
    }
	model_->setMaxLen(maxLen);

}

void PIAlign::initialize() {

    // --------------- train the base model ---------------------
    // make the model one probabilities and logify if necessary
    if(baseType_ == BASE_MODEL1 || baseType_ == BASE_MODEL1G) {
        BaseModelOne * model1 = new BaseModelOne();
        if(le2fFile_)
            model1->loadModelOne(le2fFile_,eVocab_,fVocab_,true);
        else
            model1->trainModelOne(eCorpus_, fCorpus_, eVocab_.size(), fVocab_.size(),maxSentLen_);
        if(lf2eFile_)
            model1->loadModelOne(lf2eFile_,fVocab_,eVocab_,false);
        else
            model1->trainModelOne(fCorpus_, eCorpus_, fVocab_.size(), eVocab_.size(), maxSentLen_);
        base_ = model1;
        if(baseType_ == BASE_MODEL1G)
            ((BaseModelOne*)base_)->setGeometric(true);
    }
    // train a phrase-based base measure using the dice coefficient
    else if(baseType_ == BASE_PHRASECOOC_LL) {
        base_ = new BasePhraseCooc();
        ((BasePhraseCooc*)base_)->trainCooc(eCorpus_, eVocab_, fCorpus_, fVocab_, coocDisc_, coocCut_);
    }
    // train a unigram model
    else
        base_ = new BaseUnigram();
    base_->trainUnigrams(eCorpus_, eVocab_.size(), fCorpus_, fVocab_.size());
    base_->setMaxLen(1);
    base_->trainPoisson(avgPhraseLen_, nullProb_);
    
}

// get the ID of a phrase
WordId getPhraseId(const WordString & str, StringWordMap & phrases, bool add) {
    return phrases.getId(str,add);
}
WordId getPhraseId(WordId eId, WordId fId, PairWordMap & phrases, bool add) {
    pair<WordId,WordId> myPair(eId,fId);
    return phrases.getId(myPair,add);
}

// remove a single sample
void PIAlign::removeSample(int sent) {
    // remove the sample
    model_->removeSentence(nCorpus_[sent]);
    delete nCorpus_[sent];
    nCorpus_[sent] = 0;
}

inline Prob getModelOne(const PairProbMap & model1, WordId e, WordId f) {
    PairProbMap::const_iterator it = model1.find(pair<WordId,WordId>(e,f));
    return (it==model1.end()?-50:it->second);
}

// add the generative probabilities
void PIAlign::addGenerativeProbs(const WordString & e, const WordString & f, ParseChart & chart, SpanProbMap & genChart) const {
    int maxLen = (modelType_ == MODEL_FLAT ? maxPhraseLen_ : (int)e.length());
    std::vector< LabeledEdge > eEdges = ePhrases_.findEdges(e,maxLen);
    maxLen = (modelType_ == MODEL_FLAT ? maxPhraseLen_ : (int)f.length());
    std::vector< LabeledEdge > fEdges = fPhrases_.findEdges(f,maxLen);
    for(int i = 0; i < (int)eEdges.size(); i++) {
        const LabeledEdge & ee = eEdges[i];
        for(int j = 0; j < (int)fEdges.size(); j++) {
            const LabeledEdge & fe = fEdges[j];
            PairWordMap::const_iterator it = 
                jointPhrases_.find(pair<WordId,WordId>(ee.l,fe.l));
           // cerr << " searching jointPhrases for (["<<ee.s<<","<<ee.e<<","<<ee.l<<"]/["<<fe.e<<","<<fe.s<<","<<fe.l<<"])"<<endl;
            if(it != jointPhrases_.end()) { 
                Span s(ee.s, ee.e, fe.s, fe.e);
                Prob myProb = model_->calcGenProb(it->second,s);
                if(myProb > NEG_INFINITY) {
                    chart.addToChart(s, myProb);
                    genChart.insert(SpanProbMap::value_type(s, myProb));
                }
            }
        }
    }
}

// add up the probabilities forward
void PIAlign::addForwardProbs(int eLen, int fLen, ParseChart & chart, const LookAhead & look, Prob pWidth, int hWidth, JobDetails & jd) const {

    Span yourSpan;
    Prob myProb, yourProb;
    int L = eLen+fLen,yourMax,s,t,u,v,S,U;
    if((int)jd.beamWidths.size() < L) jd.beamWidths.resize(L);
    // loop through all the agendas
    for(int l = 1; l < L; l++) {
        // get the beam and trim it to the appropriate size
        ProbSpanSet spans = chart.getTrimmedAgenda(l,hWidth,pWidth,look);
        int i, spanSize = spans.size();
        // cerr << "At length "<< l<<", processing "<< spanSize << " values"<<endl;
        for(i = 0; i < spanSize; i++) {
            const Span & mySpan = spans[i].second;
            myProb = chart.getFromChart(mySpan);
            s = mySpan.es; t = mySpan.ee; u = mySpan.fs; v = mySpan.fe;
            // cerr << "processing span "<<s<<","<<t<<","<<u<<","<<v<<endl;
            // e is forward
            for(S = max(s-l,0); S <= s; S++) {
                // f is forward
                yourMax = (S==s?u-1:u);
                for(U = max(u-l+s-S,0); U <= yourMax; U++) {
                    yourSpan = Span(S,s,U,u);
                    yourProb = chart.getFromChart(yourSpan);
                    if(yourProb > NEG_INFINITY) 
                        chart.addToChart(S,t,U,v,model_->calcTreeProb(mySpan,myProb,yourSpan,yourProb,TYPE_REG));
                }
#ifdef MONOTONIC_ON
                if(!monotonic_) {
#endif
                    // f is backward
                    yourMax = (S==s?v+1:v);
                    for(U = min(v+l-s+S, fLen); U >= yourMax; U--) {
                        yourSpan = Span(S,s,v,U);
                        yourProb = chart.getFromChart(yourSpan);
                        if(yourProb > NEG_INFINITY)
                            chart.addToChart(S,t,u,U,model_->calcTreeProb(mySpan,myProb,yourSpan,yourProb,TYPE_INV));
                    }
#ifdef MONOTONIC_ON
                }
#endif
            }
            // e is backward
            for(S = min(t+l,eLen); S >= t; S--) {
#ifdef MONOTONIC_ON
                if(!monotonic_) {
#endif
                    // f is forward
                    yourMax = (S==t?u-1:u);
                    for(U = max(u-l+S-t,0); U <= yourMax; U++) {
                        yourSpan = Span(t,S,U,u);
                        yourProb = chart.getFromChart(yourSpan);
                        if(yourProb > NEG_INFINITY)
                            chart.addToChart(s,S,U,v,model_->calcTreeProb(yourSpan,yourProb,mySpan,myProb,TYPE_INV));
                    }
#ifdef MONOTONIC_ON
                }
#endif
                // f is backward
                yourMax = (S==t?v+1:v);
                for(U = min(v+l-S+t, fLen); U >= yourMax; U--) {
                    yourSpan = Span(t,S,v,U);
                    yourProb = chart.getFromChart(yourSpan);
                    if(yourProb > NEG_INFINITY)
                        chart.addToChart(s,S,u,U,model_->calcTreeProb(yourSpan,yourProb,mySpan,myProb,TYPE_REG));
                }
            }
            // set the probability of already-covered values to zero
            //  to prevent double-adding of probabilities
            chart.removeFromChart(s,t,u,v);
        }
        // re-add the removed probabilities
        for(int j = 0; j < i; j++)
            chart[chart.findChartPosition(spans[j].second)] = spans[j].first;
        
        // cerr << "added beam of "<<i<<" at "<<l<<endl;
        jd.beamWidths[l] += i;
        jd.totalBeam += i;
        jd.totalBeamTimes++;
    }
    // std::cerr << "adding sentence" << std::endl;
    jd.sentences++;
}

void PIAlign::printSpan(const WordString & e, const WordString & f, const Span & mySpan, ostream & out, const char* phraseSep, const char* wordSep, const char* phraseBeg, const char* phraseEnd) const {
    // print if necessary
    out << phraseBeg;
    for(int i = mySpan.es; i < mySpan.ee; i++) {
        if(i != mySpan.es) out << wordSep;
        out << eVocab_.getSymbol(e[i]);
    }
    out << phraseSep;
    int sub = eVocab_.size();
    for(int i = mySpan.fs; i < mySpan.fe; i++) {
        if(i != mySpan.fs) out << wordSep;
        out << fVocab_.getSymbol(f[i]-sub);
    }
    out << phraseEnd;
}
string PIAlign::printSpan(const WordString & e, const WordString & f, const Span & mySpan, const char* phraseSep, const char* wordSep,  const char* phraseBeg, const char* phraseEnd) const {
    ostringstream oss;
    printSpan(e,f,mySpan,oss,phraseSep,wordSep,phraseBeg,phraseEnd);
    return oss.str();
}

// sample a span
SpanNode * PIAlign::sampleTree(int sent, const Span & mySpan, const ParseChart & chart, const SpanProbMap & genChart, const SpanProbMap & baseChart, bool add = true) const {
    int s=mySpan.es,t=mySpan.ee,u=mySpan.fs,v=mySpan.fe;
    // const WordString & e = eCorpus_[sent], f = fCorpus_[sent];
    // bool bracket = add;
    // reserve the probabilities
    int maxSize = (t-s+1)*(v-u+1);
    
#ifdef DEBUG_ON
    // cerr << "sampleTree("<<mySpan<<",f="<<f.length()<<",e="<<e.length()<<") == "<<printSpan(e,f,mySpan)<<endl;
    if(mySpan.length() <= 0)
        throw runtime_error("SampleTree attempted to add an empty span");
    if(maxSize <= 0)
        throw runtime_error("Maximum size less than or equal to zero");
#endif
    
    vector<Prob> probs; probs.reserve(maxSize+2);
    vector< pair<Span,Span> > pairs; pairs.reserve(maxSize);
    // Prob myProb;

    // first is the generative probability
    probs.push_back(genChart.getProb(mySpan));

    // second is the base probability
    probs.push_back(baseChart.getProb(mySpan));

    // remainder are pair probabilities
    Prob tProb;
    for(int S = s; S <= t; S++) {
        for(int U = u; U <= v; U++) {
            Span qls(s,S,u,U), qrs(S,t,U,v), rls(s,S,U,v), rrs(S,t,u,U);
            if(qls.length() && qrs.length()) {
                tProb = model_->calcTreeProb(qls,qrs,chart,TYPE_REG);
                if(tProb > NEG_INFINITY) { 
                    probs.push_back(tProb);
                    pairs.push_back(pair<Span,Span>(qls,qrs));
                }
            }
            if(rls.length() && rrs.length() && !monotonic_) {
                tProb = model_->calcTreeProb(rls,rrs,chart,TYPE_INV);
                if(tProb > NEG_INFINITY) { 
                    probs.push_back(tProb);
                    pairs.push_back(pair<Span,Span>(rls,rrs));
                }
            }
        }
    }

    // for(int i = 0; i < (int)probs.size(); i++)
    //     cerr << probs[i] << " ";
    // cerr << endl;

    // sample the answer
    int ans = sampleLogProbs(probs,annealLevel_);
    
    // make the span node
    SpanNode * myNode = new SpanNode(mySpan);
    myNode->add = add;
    myNode->prob = probs[ans];
    // if this is generative or base
    if(ans < 2) {
        // pick the type of the node
        myNode->type = ( ans == 0 ? TYPE_GEN : TYPE_BASE );
        // if not forcing word alignments, reached the bottom, or cannot proceed due to trimming, return
        if(!forceWord_ || max(mySpan.ee-mySpan.es,mySpan.fe-mySpan.fs) == 1 || probs.size() == 2)
            return myNode;
        // continue sampling
        add = false;
        ans = sampleLogProbs(&probs[2],probs.size()-2,annealLevel_)+2;
    }

    const pair<Span,Span> & myPair = pairs[ans-2];
    myNode->type = (myPair.first.fe == myPair.second.fs?TYPE_REG:TYPE_INV); 
    myNode->left = sampleTree(sent,myPair.first,chart,genChart,baseChart,add);
    myNode->right = sampleTree(sent,myPair.second,chart,genChart,baseChart,add);

    return myNode;

}

// add a single sentence sample
SpanNode * PIAlign::buildSample(int sent, ParseChart & chart, LookAhead * lookAhead, Prob pWidth, int hWidth, JobDetails & jd) const {
    timeval tStart, tInit, tBase, tGen, tLook, tFor, tSamp;
    const WordString & e = eCorpus_[sent], & f = fCorpus_[sent];
    Span sentSpan(0,e.length(),0,f.length());

#ifdef DEBUG_ON
    if(e.length() == 0 || f.length() == 0) 
        throw runtime_error("Attempted to sample zero-length sentence");
#endif
    
    // initialize
    gettimeofday(&tStart, NULL);
    chart.initialize(e.length(),f.length());
    SpanProbMap genChart = SpanProbMap();  // map of generative probs
    SpanProbMap baseChart = SpanProbMap(); // map of base probs

    // cerr << endl << "--- SAMPLING SENTENCE "<<sent<<" ---"<<endl;

    gettimeofday(&tInit, NULL);
    // add the base probabilities
    base_->addBases(e, f, *model_, chart, baseChart);
    // for(SpanProbMap::iterator it = baseChart.begin(); it != baseChart.end(); it++) {
    //     cerr << "Base: ";
    //     printSpan(e,f,it->first,cerr);
    //     cerr <<" == "<<it->second<<endl;
    // }
    gettimeofday(&tBase, NULL);

    // add the generative probabilities
    addGenerativeProbs(e,f,chart,genChart);
    gettimeofday(&tGen, NULL);

    // calculate the lookahead
    lookAhead->preCalculate(e,f,baseChart,genChart);
    gettimeofday(&tLook, NULL);

    // add the probabilities forward
    SpanNode* head;
    int eLen = e.length(), fLen = f.length();
    addForwardProbs(eLen, fLen, chart, *lookAhead, pWidth, hWidth, jd);
    Prob myLik = chart.getFromChart(sentSpan)+model_->calcSentProb(sentSpan);
    if(myLik <= NEG_INFINITY) {
        cerr << "WARNING: parsing failed! loosening beam to hist="<<histWidth_<<", prob="<<probWidth_<<endl;
        chart.setDebug(1); // base_->setDebug(1); model_->setDebug(1);
        head = buildSample(sent,chart,lookAhead,pWidth+log(10),hWidth*2,jd);
        chart.setDebug(0); // base_->setDebug(0); model_->setDebug(0);
        return head;
    }
    jd.likelihood += myLik;
    // cerr << "likelihood = "<<getFromChart(0,currELen_,0,currFLen_)<<endl;
    gettimeofday(&tFor, NULL);
    
    // sample backward probs and add sample
    head = sampleTree(sent,Span(0,eLen,0,fLen),chart,genChart,baseChart,true);

    gettimeofday(&tSamp, NULL);

    jd.timeInit += timeDifference(tStart,tInit);
    jd.timeBase += timeDifference(tInit,tBase);
    jd.timeGen += timeDifference(tBase,tGen);
    jd.timeLook += timeDifference(tGen,tLook);
    jd.timeFor += timeDifference(tLook,tFor);
    jd.timeSamp += timeDifference(tFor,tSamp);

    return head;

}

// print a single sample in tree format
void PIAlign::printSample(const WordString & e, const WordString & f, const SpanNode * myNode, ostream & sampleOut) {

    // if there are no children, print the phrase
    if(!myNode->left) 
        printSpan(e,f,myNode->span,sampleOut);
    else {
        // bracket if this is the final value generated from the actual distribution
        bool bracket = (myNode->add && myNode->left && !myNode->left->add);
        // check whether the f spans are in order or reversed
        bool ordered = (myNode->left->span.fe == myNode->right->span.fs);
        // print
        sampleOut << (bracket?"{ ":"") << (ordered?"[ ":"< ");
        printSample(e,f,myNode->left,sampleOut);
        sampleOut << " ";
        printSample(e,f,myNode->right,sampleOut);
        sampleOut << (ordered?" ]":" >") << (bracket?" }":"");
    }
        
}

inline vector<WordString> invertPhraseDic(const StringWordMap & dic) {
    std::vector<WordString> ret(dic.size(),WordString());
    for(StringWordMap::const_iterator it = dic.begin(); it != dic.end(); it++) {
        if((int)ret.size() <= it->second) ret.resize(it->second+1);
        ret[it->second] = it->first;
    }
    return ret;
}

void PIAlign::printPhraseTable(ostream & ptos) {
    vector<Prob> eProbs(ePhrases_.size(),0), fProbs(fPhrases_.size(),0), jProbs(jointPhrases_.size(),0), dProbs(jointPhrases_.size(), 0);
    model_->calcPhraseTable(jointPhrases_,eProbs,fProbs,jProbs,dProbs);
    vector<WordString> eStrs = invertPhraseDic(ePhrases_), fStrs = invertPhraseDic(fPhrases_);
    double phrasePen = exp(1);
    for(PairWordMap::const_iterator it = jointPhrases_.begin(); it != jointPhrases_.end(); it++) {
        const WordString & estr = eStrs[it->first.first];
        const WordString & fstr = fStrs[it->first.second];
        if(it->second < (int)jProbs.size() && jProbs[it->second] != 0) {
            if((int)estr.length() <= printMax_ && (int)fstr.length() <= printMax_ 
                && (int)estr.length() >= printMin_ && (int)fstr.length() >= printMin_) {
                printSpan(estr,fstr,Span(0,estr.length(),0,fstr.length()), ptos," ||| "," ","","");
                ptos << " ||| " << jProbs[it->second]/fProbs[it->first.second] <<
                        " " << jProbs[it->second]/eProbs[it->first.first] <<
                        " " << jProbs[it->second] <<
                        " " << dProbs[it->second];
                // if we are using model one, output lexical translation probabilities as well
                if(baseType_ == BASE_MODEL1 || baseType_ == BASE_MODEL1G) {
                    ptos << " " << ((BaseModelOne*)base_)->phraseModelOne(estr,fstr) << 
                            " " << ((BaseModelOne*)base_)->phraseModelOne(fstr,estr); 
                }
                ptos << " " << phrasePen << endl;
            }
        }
    }
}

template <class T>
void shuffle(vector<T> & vec) {
    int vecSize = vec.size(), pos;
    T temp;
    for(int i = vecSize - 1; i > 0; i--) {
        pos = discreteUniformSample(i+1);
        temp = vec[pos];
        vec[pos] = vec[i];
        vec[i] = temp;
    }
}

// sample all the values in the batch
void* buildSamples(void* ptr) {
    BuildJob* job = (BuildJob*)ptr;
    PIAlign * pia = job->pialign;
    for(vector<int>::iterator s = job->begin; s != job->end; s++) {
        // cerr << "Sentence "<<*s<<endl;
        SpanNode* node = pia->buildSample(*s,job->chart,job->lookAhead,pia->getProbWidth(),pia->getHistWidth(),job->details);
        pia->setNode(*s,node);
    }
    return NULL;
}

// do the whole training
void PIAlign::train() {
    
    // sample size management
    int untilNext = burnIn_, currBaby = 1, untilNextBaby = babyStepLen_,
        iter=1, currAnneal = 0, untilNextAnneal = annealStepLen_, sampNum = 0;

    // initialize parameters
    model_->sampleParameters(defStren_,defDisc_);
    vector<int> sentOrder;
    JobDetails jd;
    timeval tStart, tEnd;
    
    // iterate until we have enough samples
    while(sampNum < samples_) {

        // reset various variables
        jd.reset();
        gettimeofday(&tStart, NULL);
        // fill(derivations_.begin(), derivations_.end(), 0.0);

        // finish word iterations
        if(iter == wordIters_+1) {
            base_->setMaxLen(maxPhraseLen_);
            base_->trainPoisson(avgPhraseLen_, nullProb_);
        }

        // get info about this iteration
        if(--untilNextAnneal == 0) {
            currAnneal++;
            untilNextAnneal = annealStepLen_;
        }
        if(--untilNextBaby == 0 || sentOrder.size() == 0) {
            currBaby++;
            untilNextBaby = babyStepLen_;
            int currMaxLen = min(maxSentLen_, currBaby*maxSentLen_/babySteps_);
            // find the samples that need to be done
            for(int s = 0; s < (int)eCorpus_.size(); s++) {
                if((eCorpus_[s].length()*fCorpus_[s].length())!=0 && (int)max(eCorpus_[s].length(),fCorpus_[s].length()) <= currMaxLen) {
                    sentOrder.push_back(s);
                }
            }
        }

        // shuffle the sentence order if necessary
        if(shuffle_) shuffle(sentOrder);

        annealLevel_ = 1.0/max(1,annealSteps_-currAnneal);
        onSample_ = untilNext-- <= 0;

        cerr << "Iter "<<iter++<<" started";
        if(onSample_) cerr << ", sample "<<sampNum+1;
        cerr << endl;
        
        // main loop, sample all values
        int sents=0,lastSent=0;
        for(int beginSent = 0; beginSent < (int)sentOrder.size(); beginSent += batchLen_) {
            vector<int>::iterator beginIter = sentOrder.begin()+beginSent;
            int myBatch = min(batchLen_,(int)sentOrder.size()-beginSent);
            vector<int>::iterator endIter = beginIter+myBatch;
            // remove the samples in the current batch and count the words
            timeval tStart, tRemove;
            gettimeofday(&tStart, NULL);
            for(vector<int>::iterator s = beginIter; s != endIter; s++) {
                jd.words += eCorpus_[*s].length()+fCorpus_[*s].length();
                removeSample(*s);
                sents++;
            }
            gettimeofday(&tRemove, NULL);
            jd.timeRemove += timeDifference(tStart,tRemove);
            // cache commonly used probability values
            model_->initializeBuffers();
            // sample all the values in the batch
            for(int i = 0; i < numThreads_; i++) {
                buildJobs_[i].details.reset();
                buildJobs_[i].begin = beginIter + i*myBatch/numThreads_;
                buildJobs_[i].end = beginIter + (i+1)*myBatch/numThreads_;
                pthread_create( &buildJobs_[i].thread, NULL, buildSamples, (void*) &buildJobs_[i]);
            }
            for(int i = 0; i < numThreads_; i++) {
                pthread_join(buildJobs_[i].thread, NULL);
                jd += buildJobs_[i].details;
            }
            // for(vector<int>::iterator s = beginIter; s != endIter; s++)
            //     nCorpus_[*s] = buildSample(*s,chartTemp_);
            // add all the samples
            for(vector<int>::iterator s = beginIter; s != endIter; s++) {
                model_->addSentence(eCorpus_[*s],fCorpus_[*s],nCorpus_[*s],ePhrases_,fPhrases_,jointPhrases_);
            }
            if(sents / 100 != lastSent) {
                cerr << "\r" << sents;
                lastSent = sents/100;
            }
        }
        cerr << "\r Sentences Sampled: "<<sents<<endl;
        gettimeofday(&tEnd, NULL);
        jd.timeAll = timeDifference(tStart,tEnd);

        // print the files if we're on a sample
        if(onSample_) {
            untilNext = sampRate_-1;
            ++sampNum;
            // print the samples
            ostringstream name;
            name << prefix_ << sampNum << ".samp";
#ifdef COMPRESS_ON
            name << ".gz";
            utils::compress_ostream sampleOut(name.str().c_str(),65536);
#else
            ofstream sampleOut(name.str().c_str());
#endif
            for(int s = 0; s < (int)eCorpus_.size(); s++) {
                if(nCorpus_[s]) printSample(eCorpus_[s],fCorpus_[s],nCorpus_[s],sampleOut);
                sampleOut << endl;
            }

            // print the phrase table
            ostringstream name2;
            name2 << prefix_ << sampNum << ".pt";
#ifdef COMPRESS_ON
            name2 << ".gz";
            utils::compress_ostream ptos(name2.str().c_str(),65536);
#else
            ofstream ptos(name2.str().c_str());
#endif
            printPhraseTable(ptos);
        }

        // sample overall parameters
        model_->sampleParameters(defStren_,defDisc_);
        // tmPatterns_.sampleParameters();

        trim();

        // print stats
        model_->printStats(cerr);
        cerr << " Phrases: e="<<ePhrases_.size()<<", f="<<fPhrases_.size()<<", j="<<jointPhrases_.size()<<std::endl;
        jd.printStats(cerr);
    }

    for(int s = 0; s < (int)eCorpus_.size(); s++)
        removeSample(s);

    model_->checkEmpty();
        
} 

vector<int> phraseLengths(const StringWordMap & swm) {
    vector<int> ret;
    for(StringWordMap::const_iterator it = swm.begin(); it != swm.end(); it++) {
        if((int)ret.size() <= it->second) ret.resize(it->second+1);
        ret[it->second] = it->first.length();
    }
    return ret;
}
void PIAlign::trim() {
    vector<int> fLens = phraseLengths(fPhrases_), eLens = phraseLengths(ePhrases_);
    vector<int> fActive(fLens.size(),0), eActive(eLens.size(),0), jDead, eDead, fDead;
    // while calculating the probabilities, make sure we turn rememberNull on
    //  so we can delete only appropriate phrases
    bool remNull = model_->getRememberNull(); model_->setRememberNull(true);
    for(PairWordMap::iterator it = jointPhrases_.begin(); it != jointPhrases_.end(); it++) {
        // cerr<<"model_->calcGenProb("<<it->second<<",Span(0,"<<eLens[it->first.first]<<",0,"<<fLens[it->first.second]<<")) == "<<model_->calcGenProb(it->second,Span(0,eLens[it->first.first],0,fLens[it->first.second]))<<endl;
        if(model_->calcGenProb(it->second,Span(0,eLens[it->first.first],0,fLens[it->first.second])) > NEG_INFINITY) {
            eActive[it->first.first]++; fActive[it->first.second]++;
        }
        else
            jDead.push_back(it->second);
    }
    model_->setRememberNull(remNull);
    for(int i = 0; i < (int)eActive.size(); i++)
        if(!eActive[i] && eLens[i]) 
            eDead.push_back(i);
    for(int i = 0; i < (int)fActive.size(); i++)
        if(!fActive[i] && fLens[i]) 
            fDead.push_back(i);
    // cerr << "trimming j="<<jDead.size()<<", f="<<fDead.size()<<", e="<<eDead.size()<<endl;
    jointPhrases_.removeElements(jDead);
    ePhrases_.removeElements(eDead);
    fPhrases_.removeElements(fDead);
    
            
}
