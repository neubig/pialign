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

    void addSentence(const WordString & e, const WordString & f, SpanNode* node, StringWordMap & ePhrases, StringWordMap & fPhrases, PairWordMap & pairs);

    void removeSentence(const SpanNode* node);

};

}

#endif
