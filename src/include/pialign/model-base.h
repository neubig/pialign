#ifndef PROB_MODEL_H__
#define PROB_MODEL_H__

#include "pialign/definitions.h"
#include "pialign/pydist.h"
#include "pialign/dirichletdist.h"
#include "pialign/parse-chart.h"
#include <vector>

namespace pialign {

class ProbModel {

protected:

    // the maximum length of phrase to be generated from the base measures
    int maxLen_;

    // model probabilities
    PyDist< int,PyDenseIndex<int> > phrases_;
    DirichletDist<int> termOrNot_, forBack_;
    std::vector<int> typeCounts_;
    std::vector<Prob> derivations_;

    // buffers
    std::vector<Prob> typeProbs_;
    Prob phraseFallback_;
    int debug_;
    bool rememberNull_;


public:

    ProbModel() : phrases_(1.0,0.5), termOrNot_(1.0,2), forBack_(1.0,2),
            typeCounts_(3,0), typeProbs_(3),
            phraseFallback_(0.0), debug_(0), rememberNull_(false)  { }

    virtual ~ProbModel() { }

    // calculate the additional probability to initialize the sentence
    virtual Prob calcSentProb(const Span & mySpan) const = 0;

    // calculate the probability of generating from the generative distribution
    virtual Prob calcGenProb(WordId jId, const Span & mySpan) const = 0;
    
    // calculate the probability of generating from the base distribution
    virtual Prob calcBaseProb(const Span & mySpan, Prob baseMeas) const = 0;

    // calculate the probability of generating two trees
    virtual Prob calcTreeProb(const Span & mySpan, Prob myProb, const Span & yourSpan, Prob yourProb, int type) const = 0;
    
    // add a sentence to the distribution
    virtual Prob addSentence(const WordString & e, const WordString & f, SpanNode* node, StringWordSet & ePhrases, StringWordSet & fPhrases, PairWordSet & pairs, std::vector<Prob>& baseProbs) = 0;

    // remove a sentence sample from the distribution
    virtual SpanNode* removeSentence(const SpanNode* head, std::vector<Prob>& baseProbs) = 0;

    // initialize probabilty buffers use to save computation time
    virtual void initializeBuffers() {

        typeProbs_[0] = log(termOrNot_.getProb(0));
        typeProbs_[1] = log(termOrNot_.getProb(1)*forBack_.getProb(0));
        typeProbs_[2] = log(termOrNot_.getProb(1)*forBack_.getProb(1));
        phraseFallback_ = log(phrases_.getFallbackProb());

    }

    Prob calcTreeProb(const Span & mySpan, const Span & yourSpan, const ParseChart & chart, int type) const {
        return calcTreeProb(mySpan, chart.getFromChart(mySpan), yourSpan, chart.getFromChart(yourSpan), type);
    }

    // return the number of types that have been generated on the last turn
    //  for hierarchical models, this is zero, as type information is
    //  maintained in the tables
    const std::vector<int> & getTypeCounts() const { return typeCounts_; }

    // return the number of phrases
    virtual int getPhraseCount(WordId jId, const Span & mySpan) const {
        return phrases_.getTotal(jId);
    }

    virtual void sampleParameters(Prob defStren, Prob defDisc) {
        // std::cerr << "sampleParameters("<<defStren<<", "<<defDisc<<")"<<std::endl;
        phrases_.sampleParameters();
        if(defStren >= 0) phrases_.setStrength(defStren);
        if(defDisc >= 0) phrases_.setDiscount(defDisc);
    }
    
    virtual bool isHierarchical() = 0;

    virtual void printStats(std::ostream &out) const {
        out << " Phrase s="<<phrases_.getStrength()<<" d="<<phrases_.getDiscount()<<std::endl;
        out << " Type term="<<exp(typeProbs_[TYPE_TERM])<<" reg="<<exp(typeProbs_[TYPE_REG])<<" inv="<<exp(typeProbs_[TYPE_INV])<<std::endl;
    }

    virtual bool checkEmpty() const {
        bool ok = true;
        if(termOrNot_.getTotal() > 0) {
            std::cerr << "WARNING: non-empty termOrNot " << termOrNot_.getTotal() << std::endl;
            ok = false;
        }
        if(forBack_.getTotal() > 0) {
            std::cerr << "WARNING: non-empty forBack " << forBack_.getTotal() << std::endl;
            ok = false;
        }
        if(phrases_.getTableCount() > 0) {
           std::cerr << "WARNING: non-empty phrases " << phrases_.getTableCount() << std::endl;
           ok = false;
        }
        return ok;
    }
    
    // print the phrase table
    virtual void calcPhraseTable(const PairWordSet & jPhrases, std::vector<Prob> & eProbs, std::vector<Prob> & fProbs, std::vector<Prob> & jProbs, std::vector<Prob> & dProbs) {
        Prob myProb;
        for(PairWordSet::const_iterator it = jPhrases.begin(); it != jPhrases.end(); it++) {
            if((int)jProbs.size() <= it->second) jProbs.resize(it->second+1,0);
            myProb = phrases_.getProb(it->second,0);
            // std::cerr << "getProb("<<it->second<<") == "<<myProb<<std::endl;
            if(myProb != 0.0) {
                jProbs[it->second] = myProb;
                if((int)eProbs.size() <= it->first.first) eProbs.resize(it->first.first+1,0);
                eProbs[it->first.first] += myProb;
                if((int)fProbs.size() <= it->first.second) fProbs.resize(it->first.second+1,0);
                fProbs[it->first.second] += myProb;
            }
        }
        dProbs = derivations_;
    }

    Prob addAverageDerivation(WordId jId, int count, Prob prob) {
#ifdef DEBUG_ON
        if(prob == 0.0) THROW_ERROR("Attempting to add 0-valued derivation @ "<<jId);
#endif
        if((int)derivations_.size() <= jId) derivations_.resize(jId+1,0);
        // std::cerr << "derivations_["<<jId<<"] = ("<<derivations_[jId]<<"*("<<count-1<<") + exp("<<prob<<"))/"<<count<<std::endl;
        derivations_[jId] = (derivations_[jId]*(count-1) + exp(prob))/count;
        return derivations_[jId];
    }

    // add a type and return its log probability
    inline Prob addType(WordId id) {
        PRINT_DEBUG(" addType("<<id<<")");
        Prob ret = 0;
        if(id == 0) {
            ret = log(termOrNot_.getProb(0));
            termOrNot_.add(0);
        } else {
            ret = log(termOrNot_.getProb(1)*forBack_.getProb(id-1));
            termOrNot_.add(1); forBack_.add(id-1);
        }
        PRINT_DEBUG(" == "<<ret<<std::endl);
        return ret;
    }

    // remove a type and return its log probability
    inline Prob removeType(WordId id) {
        PRINT_DEBUG(" removeType("<<id<<")");
        Prob ret;
        if(id == 0) {
            termOrNot_.remove(0);
            ret = log(termOrNot_.getProb(0));
        } else {
            termOrNot_.remove(1); forBack_.remove(id-1);
            ret = log(termOrNot_.getProb(1)*forBack_.getProb(id-1));
        }
        PRINT_DEBUG(" == "<<ret<<std::endl);
        return ret;
    }

	virtual void setMaxLen(int maxLen) { maxLen_ = maxLen; }

    bool getRememberNull() { return rememberNull_; }
    void setRememberNull(bool rememberNull) { rememberNull_ = rememberNull; }
    void setDebug(int debug) { debug_ = debug; }
    void setTerminalPrior(Prob termPrior) {
        std::vector<Prob> vec(2); vec[0] = termPrior; vec[1] = 1-termPrior;
        termOrNot_.setBases(vec);
    }
    void setTerminalStrength(Prob stren) { termOrNot_.setStrength(stren); }

};

}

#endif
