#include "pialign.h"

// #include "statecollection.h"
#include <fstream>
#include <cmath>
#include <sys/time.h>


#include "base-model1.h"
#include "base-unigram.h"

#include "prob-hier.h"
#include "prob-flat.h"
#include "prob-length.h"

#ifdef COMPRESS_ON
#include "compress_stream.hpp"
#endif

#ifdef HAVE_CONFIG_H 
#include "config.h"
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
<< "               'uni'=simple unigrams (default 'm1g')" << endl
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
<< " -samps        The number of samples to take (default 1)" << endl
<< " -samprate     Take samples every samprate turns (default 1)" << endl
<< " -worditers    The number of iterations to perform with a word-based model (default 0)" << endl << endl;
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
            else if(!strcmp(argv[i],"-noword"))         forceWord_ = false;
            else if(!strcmp(argv[i],"-noremnull"))      rememberNull = false;
            else if(!strcmp(argv[i],"-monotonic"))      monotonic_ = true;
            else if(!strcmp(argv[i],"-babysteps"))      babySteps_ = atoi(argv[++i]);
            else if(!strcmp(argv[i],"-babysteplen"))    babyStepLen_ = atoi(argv[++i]);
            else if(!strcmp(argv[i],"-annealsteps"))    annealSteps_ = atoi(argv[++i]);
            else if(!strcmp(argv[i],"-annealsteplen"))  annealStepLen_ = atoi(argv[++i]);
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
                else {
                    ostringstream oss;
                    oss << "Unknown base argument "<<argv[i];
                    dieOnHelp(oss.str().c_str());
                }
            }
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

}

// load the corpus
Corpus PIAlign::loadCorpus(string file, SymbolSet< string, WordId > & vocab, WordId boost) {
    ifstream ifs(file.c_str());
    if(!ifs) {
        ostringstream oss;
        oss << "IO Error: Couldn't open "<<file;
        throw runtime_error(oss.str());
    }
    string line,buff;
    Corpus ret;
    while(getline(ifs, line)) {
        istringstream iss(line);
        vector<WordId> sent;
        while(iss >> buff)
            sent.push_back(vocab.getId(buff,true)+boost);
        ret.push_back(GenericString<WordId>(sent));
    }
    return ret;
}

void PIAlign::loadCorpora() {

    // load the corpora
    eCorpus_ = loadCorpus(eFile_, eVocab_, 0);
    fCorpus_ = loadCorpus(fFile_, fVocab_, eVocab_.size());
    // remove sentences that are too long
    int total = 0, maxLen = 0;
    for(unsigned i = 0; i < eCorpus_.size(); i++) {
        int myLen = max(eCorpus_[i].length(),fCorpus_[i].length());
        if(myLen > maxSentLen_) {
            eCorpus_[i] = WordString();
            fCorpus_[i] = WordString();
        } else {
            total++;
            maxLen = max(myLen,maxLen);
        }
    }
    if(eCorpus_.size() != fCorpus_.size())
       throw std::runtime_error("Corpus sizes don't match");
    cerr << "Loaded corpus: "<<total<<"/"<<eCorpus_.size()<<" sentences used"<<endl;
    
    // allocate other corpora
    tCorpus_ = vector<int>(eCorpus_.size()*3,0);
    hCorpus_ = vector<int>(eCorpus_.size(),-1);
    aCorpus_ = Corpus(eCorpus_.size(), WordString());
    
    // initialize the chart
    chart_.initialize(maxLen,maxLen);
	model_->setMaxLen(maxLen);

}

void PIAlign::initialize() {

    // make the model one probabilities and logify if necessary
    if(baseType_ != BASE_UNI) {
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
    } else
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
    model_->removeSentence(hCorpus_[sent], aCorpus_[sent], &tCorpus_[sent*3]);
}

inline Prob getModelOne(const PairProbMap & model1, WordId e, WordId f) {
    PairProbMap::const_iterator it = model1.find(pair<WordId,WordId>(e,f));
    return (it==model1.end()?-50:it->second);
}

// find all active phrases in a string
// TODO: change this to use Tries
vector<LabeledEdge> PIAlign::findEdges(const WordString & str, const StringWordMap & dict) {
    int T = str.length();
    int maxLen = (modelType_==MODEL_FLAT?maxPhraseLen_:T);
    vector<LabeledEdge> ret;    
    for(int i = 0; i <= T; i++) {
        for(int j = i; j <= min(i+maxLen,T); j++) {
            StringWordMap::const_iterator it = dict.find(str.substr(i,j-i));
            if(it != dict.end()) {
                ret.push_back(LabeledEdge(i,j,it->second));
            }
        }
    }
    return ret;
}

// add the generative probabilities
void PIAlign::addGenerativeProbs(const WordString & e, const WordString & f) {
    std::vector< LabeledEdge > eEdges = findEdges(e,ePhrases_);
    std::vector< LabeledEdge > fEdges = findEdges(f,fPhrases_);
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
                    chart_.addToChart(s, myProb);
                    genChart_.insert(SpanProbMap::value_type(s, myProb));
                }
            }
        }
    }
}

// add up the probabilities forward
void PIAlign::addForwardProbs(int eLen, int fLen) {

    Span yourSpan;
    Prob myProb, yourProb;
    int L = eLen+fLen,yourMax,s,t,u,v,S,U;
    // loop through all the agendas
    for(int l = 1; l < L; l++) {
        // get the beam and trim it to the appropriate size
        ProbSpanSet spans = chart_.getTrimmedAgenda(l,histWidth_,probWidth_);
        int i, spanSize = spans.size();
        // cerr << "At length "<< l<<", processing "<< spanSize << " values"<<endl;
        for(i = 0; i < spanSize; i++) {
            const Span & mySpan = spans[i].second;
            myProb = spans[i].first;
            s = mySpan.es; t = mySpan.ee; u = mySpan.fs; v = mySpan.fe;
            // cerr << "processing span "<<s<<","<<t<<","<<u<<","<<v<<endl;
            // e is forward
            for(S = max(s-l,0); S <= s; S++) {
                // f is forward
                yourMax = (S==s?u-1:u);
                for(U = max(u-l+s-S,0); U <= yourMax; U++) {
                    yourSpan = Span(S,s,U,u);
                    yourProb = chart_.getFromChart(yourSpan);
                    if(yourProb > NEG_INFINITY) 
                        chart_.addToChart(S,t,U,v,model_->calcTreeProb(mySpan,myProb,yourSpan,yourProb,TYPE_REG));
                }
#ifdef MONOTONIC_ON
                if(!monotonic_) {
#endif
                    // f is backward
                    yourMax = (S==s?v+1:v);
                    for(U = min(v+l-s+S, fLen); U >= yourMax; U--) {
                        yourSpan = Span(S,s,v,U);
                        yourProb = chart_.getFromChart(yourSpan);
                        if(yourProb > NEG_INFINITY)
                            chart_.addToChart(S,t,u,U,model_->calcTreeProb(mySpan,myProb,yourSpan,yourProb,TYPE_INV));
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
                        yourProb = chart_.getFromChart(yourSpan);
                        if(yourProb > NEG_INFINITY)
                            chart_.addToChart(s,S,U,v,model_->calcTreeProb(yourSpan,yourProb,mySpan,myProb,TYPE_INV));
                    }
#ifdef MONOTONIC_ON
                }
#endif
                // f is backward
                yourMax = (S==t?v+1:v);
                for(U = min(v+l-S+t, fLen); U >= yourMax; U--) {
                    yourSpan = Span(t,S,v,U);
                    yourProb = chart_.getFromChart(yourSpan);
                    if(yourProb > NEG_INFINITY)
                        chart_.addToChart(s,S,u,U,model_->calcTreeProb(yourSpan,yourProb,mySpan,myProb,TYPE_REG));
                }
            }
            // set the probability of already-covered values to zero
            //  to prevent double-adding of probabilities
            chart_.removeFromChart(s,t,u,v);
        }
        // re-add the removed probabilities
        for(int j = 0; j < i; j++)
            chart_[chart_.findChartPosition(spans[j].second)] = spans[j].first;
        totalBeam_ += i;
        totalBeamTimes_++;
    }
}

void PIAlign::printSpan(const WordString & e, const WordString & f, const Span & mySpan, ostream & out, const char* phraseSep, const char* wordSep, const char* phraseBeg, const char* phraseEnd) {
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
string PIAlign::printSpan(const WordString & e, const WordString & f, const Span & mySpan, const char* phraseSep, const char* wordSep,  const char* phraseBeg, const char* phraseEnd) {
    ostringstream oss;
    printSpan(e,f,mySpan,oss,phraseSep,wordSep,phraseBeg,phraseEnd);
    return oss.str();
}

Prob PIAlign::getGenProb(const Span & mySpan) {
    SpanProbMap::const_iterator it = genChart_.find(mySpan);
    return it == genChart_.end() ? NEG_INFINITY : it->second;
}

// sample a span
WordId PIAlign::sampleTree(int sent, const Span & mySpan, WordString & sentIds, bool add = true) {
    int s=mySpan.es,t=mySpan.ee,u=mySpan.fs,v=mySpan.fe;
    const WordString & e = eCorpus_[sent], f = fCorpus_[sent];
    bool bracket = add;
#ifdef DEBUG_ON
    // cerr << "sampleTree("<<mySpan<<",f="<<f.length()<<",e="<<e.length()<<") == "<<printSpan(e,f,mySpan)<<endl;
    if(mySpan.length() == 0)
        throw runtime_error("SampleTree attempted to add an empty span");
#endif

    // reserve the probabilities
    int maxSize = (t-s+1)*(v-u+1);
    vector<Prob> probs; probs.reserve(maxSize+2);
    vector< pair<Span,Span> > pairs; pairs.reserve(maxSize);
    // Prob myProb;

    // first is the generative probability
    probs.push_back(getGenProb(mySpan));

    // second is the base probability
    probs.push_back(base_->getFromChart(mySpan));

    // remainder are pair probabilities
    Prob tProb;
    for(int S = s; S <= t; S++) {
        for(int U = u; U <= v; U++) {
            Span qls(s,S,u,U), qrs(S,t,U,v), rls(s,S,U,v), rrs(S,t,u,U);
            if(qls.length() && qrs.length()) {
                tProb = model_->calcTreeProb(qls,qrs,chart_,TYPE_REG);
                if(tProb > NEG_INFINITY) { 
                    probs.push_back(tProb);
                    pairs.push_back(pair<Span,Span>(qls,qrs));
                }
            }
            if(rls.length() && rrs.length() && !monotonic_) {
                tProb = model_->calcTreeProb(rls,rrs,chart_,TYPE_INV);
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
    
    // generative or base dist
    WordId eId = getPhraseId(eCorpus_[sent].substr(s,t-s),ePhrases_,add),
        fId = getPhraseId(fCorpus_[sent].substr(u,v-u),fPhrases_,add),
        jId = -1;
    if(ans < 2) {
        jId = getPhraseId(eId,fId,jointPhrases_,add);
        if(add) {
            if(ans == 0)
                model_->addGen(jId,mySpan,probs[0]);
            else
                model_->addBase(jId, mySpan,probs[1]);
        }
        if(!sampleOut_ || !forceWord_ || max(mySpan.ee-mySpan.es,mySpan.fe-mySpan.fs) == 1) {
            if(sampleOut_) printSpan(e,f,mySpan,*sampleOut_);
            return jId;
        }
        add = false;
        ans = sampleLogProbs(&probs[2],probs.size()-2,annealLevel_)+2;
    }
    else
        bracket = false;

    const pair<Span,Span> & myPair = pairs[ans-2];
    int type = (myPair.first.fe == myPair.second.fs?TYPE_REG:TYPE_INV);
    if(sampleOut_) *sampleOut_ << (bracket?"{ ":"") << (type==TYPE_REG?"[ ":"< ");
    WordId lId = sampleTree(sent,myPair.first,sentIds,add);
    if(sampleOut_) *sampleOut_ << " ";
    WordId rId = sampleTree(sent,myPair.second,sentIds,add);
    if(sampleOut_) *sampleOut_ << (type==TYPE_REG?" ]":" >") << (bracket?" }":"");
    
    if(add && model_->isHierarchical()) {
        jId = getPhraseId(eId,fId,jointPhrases_,add);
        model_->addTree(jId,lId,rId,mySpan,myPair.first,myPair.second,type,probs[ans]);
    }

    return jId;

}

double timeDifference(const timeval & s, const timeval & e) {
    return (e.tv_sec-s.tv_sec)+(e.tv_usec-s.tv_usec)/1000000.0;
}

// add a single sentence sample
void PIAlign::addSample(int sent) {
    timeval tStart, tRemove, tInit, tBase, tGen, tFor, tSamp;
    const WordString & e = eCorpus_[sent], & f = fCorpus_[sent];
    Span sentSpan(0,e.length(),0,f.length());

    // remove the sample
    gettimeofday(&tStart, NULL);
    removeSample(sent);
    genChart_ = SpanProbMap();
    
    // initialize
    gettimeofday(&tRemove, NULL);
    chart_.initialize(e.length(),f.length());
    base_->initialize(e.length(),f.length());
    model_->initialize(e,f,chart_,aCorpus_[sent],&tCorpus_[sent*3]);

    // cerr << endl << "--- SAMPLING SENTENCE "<<sent<<" ---"<<endl;

    gettimeofday(&tInit, NULL);
    // add the base probabilities
    base_->addBases(e, f, *model_, chart_);
    gettimeofday(&tBase, NULL);

    // add the generative probabilities
    addGenerativeProbs(e,f);
    gettimeofday(&tGen, NULL);

    // add the probabilities forward
    int eLen = e.length(), fLen = f.length();
    addForwardProbs(eLen, fLen);
    Prob myLik = chart_.getFromChart(sentSpan)+model_->calcSentProb(sentSpan);
    if(myLik <= NEG_INFINITY) {
        probWidth_ += log(10); histWidth_ *= 2;
        cerr << "WARNING: parsing failed! loosening beam to hist="<<histWidth_<<", prob="<<probWidth_<<endl;
        chart_.setDebug(1); base_->setDebug(1); model_->setDebug(1);
        addSample(sent);
        chart_.setDebug(0); base_->setDebug(0); model_->setDebug(0);
        probWidth_ -= log(10); histWidth_ /= 2;
        return;
    }
    likelihood_ += myLik;
    // cerr << "likelihood = "<<getFromChart(0,currELen_,0,currFLen_)<<endl;
    gettimeofday(&tFor, NULL);
    
    // sample backward probs and add sample
    hCorpus_[sent] = sampleTree(sent,Span(0,eLen,0,fLen),aCorpus_[sent],true);
    if(sampleOut_) *sampleOut_ << endl;

    gettimeofday(&tSamp, NULL);

    timeRemove_ += timeDifference(tStart,tRemove);
    timeInit_ += timeDifference(tRemove,tInit);
    timeBase_ += timeDifference(tInit,tBase);
    timeGen_ += timeDifference(tBase,tGen);
    timeFor_ += timeDifference(tGen,tFor);
    timeSamp_ += timeDifference(tFor,tSamp);
    timeAll_ += timeDifference(tStart,tSamp);

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
                if(baseType_ != BASE_UNI) {
                    ptos << " " << ((BaseModelOne*)base_)->phraseModelOne(estr,fstr) << 
                            " " << ((BaseModelOne*)base_)->phraseModelOne(fstr,estr); 
                }
                ptos << " " << phrasePen << endl;
            }
        }
    }
}

// do the whole training
void PIAlign::train() {
    
    // sample size management
    int untilNext = burnIn_, currBaby = 1, untilNextBaby = babyStepLen_,
        iter=1, currAnneal = 0, untilNextAnneal = annealStepLen_, sampNum = 0;

    // initialize parameters
    model_->sampleParameters(defStren_,defDisc_);
    
    // iterate until we have enough samples
    while(sampNum < samples_) {

        // set various variables to zero
        likelihood_ = 0;
        timeRemove_ = 0; timeInit_ = 0; timeBase_ = 0; timeGen_ = 0;
        timeFor_ = 0; timeSamp_ = 0; timeAll_ = 0;
        totalBeam_ = 0, totalBeamTimes_ = 0;
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
        if(--untilNextBaby == 0) {
            currBaby++;
            untilNextBaby = babyStepLen_;
        }
        int currMaxLen = min(maxSentLen_, currBaby*maxSentLen_/babySteps_);
        annealLevel_ = 1.0/max(1,annealSteps_-currAnneal);
        onSample_ = untilNext-- <= 0;
        sampleOut_ = 0;
        if(onSample_) {
            untilNext = sampRate_-1;
            ++sampNum;
            ostringstream name;
            name << prefix_ << sampNum << ".samp";
#ifdef COMPRESS_ON
            name << ".gz";
            sampleOut_ = new utils::compress_ostream(name.str().c_str(),65536);
#else
            sampleOut_ = new ofstream(name.str().c_str());
#endif
        }
        likelihood_ = 0;
       cerr << "Iter "<<iter++<<" started: a="<<annealLevel_<<", ms="<<currMaxLen<<", mp="<<base_->getMaxLen();
        if(onSample_) cerr << ", sample "<<sampNum;
       cerr << endl;
        
        // sample values, break when necessary
        int words=0,sents=0;
        for(int s = 0; s < (int)eCorpus_.size(); s++) {
            if((eCorpus_[s].length()*fCorpus_[s].length())!=0 && (int)max(eCorpus_[s].length(),fCorpus_[s].length()) <= currMaxLen) {
                words += eCorpus_[s].length()+fCorpus_[s].length();
                addSample(s);
                if(++sents % 100 == 0) cerr << "\r" << sents;
            }
            else if(onSample_)
                *sampleOut_ << endl;
        }
       cerr << "\r Sentences Sampled: "<<sents<<endl;

        if(onSample_) {
            delete sampleOut_;
            sampleOut_ = 0;
            ostringstream name;
            name << prefix_ << sampNum << ".pt";
#ifdef COMPRESS_ON
            name << ".gz";
            utils::compress_ostream ptos(name.str().c_str(),65536);
#else
            ofstream ptos(name.str().c_str());
#endif
            printPhraseTable(ptos);
        }

        // sample overall parameters
        model_->sampleParameters(defStren_,defDisc_);
        // tmPatterns_.sampleParameters();

        trim();

        // print stats
        model_->printStats(cerr);
        cerr << " Likelihood="<< likelihood_/words <<endl;       
        cerr << " Phrases: e="<<ePhrases_.size()<<", f="<<fPhrases_.size()<<", j="<<jointPhrases_.size()<<endl;
        cerr << " Time="<<timeAll_<<"s (r="<<timeRemove_<<", i="<<timeInit_<<", b="<<timeBase_<<", g="<<timeGen_<<", f="<<timeFor_<<", s="<<timeSamp_<<")"<<endl;
        cerr << " Avg. Beam="<<(double)totalBeam_/totalBeamTimes_<<endl;

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
