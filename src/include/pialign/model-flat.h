#ifndef PROB_FLAT_H__
#define PROB_FLAT_H__

#include "model-base.h"

namespace pialign {

class FlatModel : public ProbModel {
    
public:

    FlatModel() : ProbModel()
     { }

    bool isHierarchical() { return false; }
    
    Prob calcSentProb(const Span & mySpan) const { return 0; }

    inline bool isNull(const Span & mySpan) const {
        return (mySpan.ee == mySpan.es || mySpan.fe == mySpan.fs);
    }

    Prob calcGenProb(WordId jId, const Span & mySpan) const {
        return typeProbs_[TYPE_TERM]+log(phrases_.getProb(jId,0));
    }
    
    Prob calcBaseProb(const Span & mySpan, Prob baseMeas) const {
#ifdef DEBUG_ON
        if(debug_)
            std::cerr << "FlatModel::calcBaseProb @ "<<mySpan<<" ("<<typeProbs_[TYPE_TERM]<<"+"<<((!isNull(mySpan)||rememberNull_)?phraseFallback_:0)<<"+"<<baseMeas<<") == "<<typeProbs_[TYPE_TERM]+baseMeas+((!isNull(mySpan)||rememberNull_)?phraseFallback_:0)<<std::endl;
#endif
        return typeProbs_[TYPE_TERM]+baseMeas+(!isNull(mySpan)||rememberNull_?phraseFallback_:0);
    }
    Prob calcTreeProb(const Span & mySpan, Prob myProb, const Span & yourSpan, Prob yourProb, int type) const {
        return typeProbs_[type]+myProb+yourProb;
    }
    void addGen(WordId jId, const Span & mySpan, Prob prob, int* tCounts, WordString & jIds) {
        tCounts[TYPE_TERM]++;
        addType(TYPE_TERM);
        if(rememberNull_ || !isNull(mySpan)) {
            jIds = jIds + jId;
            phrases_.addExisting(jId);
            addAverageDerivation(jId,phrases_.getTotal(jId),prob);
        } 
    }

    void addBase(WordId jId, const Span & mySpan, Prob prob, int* tCounts, WordString & jIds) {
        tCounts[TYPE_TERM]++;
        addType(TYPE_TERM);
        if(rememberNull_ || !isNull(mySpan)) {
            jIds = jIds + jId;
            phrases_.addNew(jId,-1,-1,-1);
            addAverageDerivation(jId,phrases_.getTotal(jId),prob);
        }
    }

    void addTree(WordId jId, WordId lId, WordId rId, 
                const Span & jSpan, const Span & lSpan, const Span & rSpan, 
                int type, Prob prob, int* tCounts, WordString & jIds) {
        tCounts[type]++;
        addType(type);
    }
    void removeSentence(WordId head, WordString & jIds, int* counts) {
        for(int i = 0; i < (int)jIds.length(); i++) phrases_.remove(jIds[i]);
        jIds = WordString();
        for(int i = 0; i < 3; i++) 
            for( ; counts[i] != 0; counts[i]--) 
                removeType(i);
    }

};

}

#endif
