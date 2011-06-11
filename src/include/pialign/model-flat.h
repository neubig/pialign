#ifndef PROB_FLAT_H__
#define PROB_FLAT_H__

#include "pialign/model-base.h"

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
        PRINT_DEBUG("FlatModel::calcBaseProb @ "<<mySpan<<" ("<<typeProbs_[TYPE_TERM]<<"+"<<((!isNull(mySpan)||rememberNull_)?phraseFallback_:0)<<"+"<<baseMeas<<") == "<<typeProbs_[TYPE_TERM]+baseMeas+((!isNull(mySpan)||rememberNull_)?phraseFallback_:0)<<std::endl);
        return typeProbs_[TYPE_TERM]+baseMeas+(!isNull(mySpan)||rememberNull_?phraseFallback_:0);
    }
    Prob calcTreeProb(const Span & mySpan, Prob myProb, const Span & yourSpan, Prob yourProb, int type) const {
        PRINT_DEBUG("FlatModel::calcTree("<<mySpan<<"/"<<yourSpan<<") == " <<typeProbs_[type]<<"+"<<myProb<<"+"<<yourProb<<" == "<<typeProbs_[type]+myProb+yourProb<<std::endl);
        return typeProbs_[type]+myProb+yourProb;
    }

    Prob addSentence(const WordString & e, const WordString & f, SpanNode* node, StringWordSet & ePhrases, StringWordSet & fPhrases, PairWordSet & pairs, std::vector<Prob>& baseProbs);

    SpanNode* removeSentence(const SpanNode* node, std::vector<Prob>& baseProbs);

};

}

#endif
