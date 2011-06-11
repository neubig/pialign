#ifndef PROB_LENGTH_H__
#define PROB_LENGTH_H__

#include "pialign/model-base.h"
#include <stdexcept>

namespace pialign {

class LengthModel : public ProbModel {

protected:

    std::vector< PyDist< WordId,PySparseIndex<WordId> > > sepPhrases_;
    std::vector< DirichletDist<int> > sepFor_, sepTerm_;
    std::vector< std::vector<Prob> > sepType_;
    std::vector< Prob > sepFallbacks_, sepSplits_;
    std::vector< int > phraseIdxs_, hist;
    Prob sentPen_;

public:

    LengthModel() : ProbModel() {
		
	}
	
	void setMaxLen(int maxLen);
    
    Prob calcSentProb(const Span & mySpan) const { return sentPen_; }
    bool isHierarchical() { return true; }

    Prob calcGenProb(WordId jId, const Span & mySpan) const {
        int idx = mySpan.length()-1;
        return (idx || rememberNull_) ? log(sepPhrases_[idx].getProb(jId,0)) : NEG_INFINITY;
    }
    
    Prob calcBaseProb(const Span & mySpan, Prob baseMeas) const {
        int idx = mySpan.length()-1;
        PRINT_DEBUG("LengthModel::calcBaseProb"<<mySpan<<" == "<<sepFallbacks_[idx]<<"+"<<sepType_[idx][TYPE_TERM]<<"+"<<baseMeas<<" == "<<sepFallbacks_[idx]+sepType_[idx][TYPE_TERM]+baseMeas<<std::endl);
        return sepFallbacks_[idx]+sepType_[idx][TYPE_TERM]+baseMeas;
    }

    Prob calcTreeProb(const Span & mySpan, Prob myProb, const Span & yourSpan, Prob yourProb, int type) const {
        int idx = mySpan.length()+yourSpan.length()-1;
#ifdef DEBUG_ON
        if(idx == 0)
            throw std::runtime_error("illegal value in calcTreeProb");
#endif
        PRINT_DEBUG("LengthModel::calcTreeProb("<<mySpan<<","<<myProb<<","<<yourSpan<<","<<yourProb<<","<<type<<") == "<<sepFallbacks_[idx]<<"+"<<sepSplits_[idx]<<"+"<<sepType_[idx][type]<<"+"<<myProb<<"+"<<yourProb<<" == "<<sepFallbacks_[idx]+sepSplits_[idx]+sepType_[idx][type]+myProb+yourProb<<std::endl);
        return sepFallbacks_[idx]+sepSplits_[idx]+sepType_[idx][type]+myProb+yourProb;
    }
    inline int saveIdx(WordId jId, int index) {
        while((int)phraseIdxs_.size() <= jId) 
            phraseIdxs_.resize(jId+1);
        return (phraseIdxs_[jId] = index);
    }

    
    Prob addSentence(const WordString & e, const WordString & f, SpanNode* node, StringWordSet & ePhrases, StringWordSet & fPhrases, PairWordSet & pairs, std::vector<Prob>& baseProbs);

    SpanNode* removePhrasePair(WordId jId, std::vector<Prob>& baseProbs);

    SpanNode* removeSentence(const SpanNode* node, std::vector<Prob>& baseProbs) {
        if(!node) return 0;
        return removePhrasePair(node->phraseid, baseProbs);
    }
    
    void initialize(const WordString & e, const WordString & f, 
            ParseChart & chart, WordString & jIds, int* tCounts);
    
    void sampleParameters(Prob defStren, Prob defDisc) {
        for(int i = 0; i < (int)sepPhrases_.size(); i++) {
            sepPhrases_[i].sampleParameters();
            if(defStren >= 0) sepPhrases_[i].setStrength(defStren);
            if(defDisc >= 0) sepPhrases_[i].setDiscount(defDisc);
        }
    }

    virtual void printStats(std::ostream &out) const;

    void calcPhraseTable(const PairWordSet & jPhrases, std::vector<Prob> & eProbs, std::vector<Prob> & fProbs, std::vector<Prob> & jProbs, std::vector<Prob> & dProbs);

    inline Prob addType(WordId id, int idx) {
        Prob ret;
        if(id == 0) {
            ret = log(sepTerm_[idx].getProb(0));
            sepTerm_[idx].add(0);
        }
        else {
            ret = log(sepTerm_[idx].getProb(1)*sepFor_[idx].getProb(id-1));
            sepTerm_[idx].add(1); sepFor_[idx].add(id-1);
        }
        // std::cerr << "addType("<<id<<","<<idx<<") == "<<ret<<std::endl;
        return ret;
    }
    inline Prob removeType(WordId id, int idx) {
        Prob ret;
        if(id == 0) {
            sepTerm_[idx].remove(0);
            ret = log(sepTerm_[idx].getProb(0));
        } else {
            sepTerm_[idx].remove(1); sepFor_[idx].remove(id-1);
            ret = log(sepTerm_[idx].getProb(1)*sepFor_[idx].getProb(id-1));
        }
        // std::cerr << "removeType("<<id<<","<<idx<<") == "<<ret<<std::endl;
        return ret;
    }

    // initialize probabilty buffers use to save computation time
    void initializeBuffers() {

        for(unsigned idx = 0; idx < sepType_.size(); idx++) {
            sepType_[idx][0] = log(sepTerm_[idx].getProb(0));
            sepType_[idx][1] = log(sepTerm_[idx].getProb(1)*sepFor_[idx].getProb(0));
            sepType_[idx][2] = log(sepTerm_[idx].getProb(1)*sepFor_[idx].getProb(1));
            sepFallbacks_[idx] = log(sepPhrases_[idx].getFallbackProb());
        }

    }

};

}

#endif
