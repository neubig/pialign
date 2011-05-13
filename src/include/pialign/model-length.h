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
        return sepFallbacks_[idx]+sepType_[idx][TYPE_TERM]+baseMeas;
    }

    Prob calcTreeProb(const Span & mySpan, Prob myProb, const Span & yourSpan, Prob yourProb, int type) const {
        int idx = mySpan.length()+yourSpan.length()-1;
#ifdef DEBUG_ON
        if(idx == 0)
            throw std::runtime_error("illegal value in calcTreeProb");
#endif
        return sepFallbacks_[idx]+sepSplits_[idx]+sepType_[idx][type]+myProb+yourProb;
    }
    inline int saveIdx(WordId jId, int index) {
        while((int)phraseIdxs_.size() <= jId) 
            phraseIdxs_.resize(jId+1);
        return (phraseIdxs_[jId] = index);
    }

    // void addGen(WordId jId, const Span & mySpan, Prob prob) {
    //     int idx = saveIdx(jId,mySpan.length()-1);
    //     // std::cerr << "addGen("<<jId<<") --> "<<idx<<std::endl;
    //     sepPhrases_[idx].addExisting(jId);
    //     addAverageDerivation(jId,sepPhrases_[idx].getTotal(jId),prob);
    // }

    // void addBase(WordId jId, const Span & mySpan, Prob prob) {
    //     int idx = saveIdx(jId,mySpan.length()-1);
    //     // std::cerr << "addBase("<<jId<<") --> "<<idx<<std::endl;
    //     sepPhrases_[idx].addNew(jId,-1,-1,TYPE_TERM);
    //     addAverageDerivation(jId,sepPhrases_[idx].getTotal(jId),prob);
    //     addType(TYPE_TERM,idx);
    // }

    // void addTree(WordId jId, WordId lId, WordId rId, 
    //             const Span & jSpan, const Span & lSpan, const Span & rSpan, 
    //             int type, Prob prob) {
    //     int idx = saveIdx(jId,jSpan.length()-1);
    //     // std::cerr << "addTree("<<jId<<") --> "<<idx<<std::endl;
    //     sepPhrases_[idx].addNew(jId, lId, rId, type);
    //     addAverageDerivation(jId,sepPhrases_[idx].getTotal(jId),prob);
    //     addType(type,idx);
    // }
    
    void addSentence(const WordString & e, const WordString & f, SpanNode* node, StringWordMap & ePhrases, StringWordMap & fPhrases, PairWordMap & pairs);

    void removePhrasePair(WordId jId);

    void removeSentence(const SpanNode* node) {
        if(!node) return;
        removePhrasePair(node->phraseid);
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

    void calcPhraseTable(const PairWordMap & jPhrases, std::vector<Prob> & eProbs, std::vector<Prob> & fProbs, std::vector<Prob> & jProbs, std::vector<Prob> & dProbs);

    inline void addType(WordId id, int idx) {
        if(id == 0)
            sepTerm_[idx].add(0);
        else {
            sepTerm_[idx].add(1);
            sepFor_[idx].add(id-1);
        }
    }
    inline void removeType(WordId id, int idx) {
        if(id == 0)
            sepTerm_[idx].remove(0);
        else {
            sepTerm_[idx].remove(1);
            sepFor_[idx].remove(id-1);
        }
    }

};

}

#endif
