#ifndef PROB_HIER_H__
#define PROB_HIER_H__

#include "pialign/model-base.h"
#include <stdexcept>

namespace pialign {

class HierModel : public ProbModel {

public:

    HierModel() : ProbModel() { }
    
    Prob calcSentProb(const Span & mySpan) const { return 0; }

    Prob calcGenProb(WordId jId, const Span & mySpan) const {
        return log(phrases_.getProb(jId,0));
    }
    
    Prob calcBaseProb(const Span & mySpan, Prob baseMeas) const {
#ifdef DEBUG_ON
        if(debug_)
            std::cerr << "HierModel::calcBaseProb("<<phraseFallback_<<"+"<<typeProbs_[TYPE_TERM]<<"+"<<baseMeas<<")"<<std::endl;
#endif
        return phraseFallback_+typeProbs_[TYPE_TERM]+baseMeas;
    }
    
    Prob calcTreeProb(const Span & mySpan, Prob myProb, const Span & yourSpan, Prob yourProb, int type) const {
        // std::cerr << " HierModel::calcTreeProb @ "<<mySpan<<"/"<<yourSpan<<" == ("<<phraseFallback_<<","<<typeProbs_[type]<<","<<myProb<<","<<yourProb<<")"<<std::endl;
        return phraseFallback_+typeProbs_[type]+myProb+yourProb;
    }

    void addSentence(const WordString & e, const WordString & f, SpanNode* node, StringWordMap & ePhrases, StringWordMap & fPhrases, PairWordMap & pairs);

    void removeSentence(const SpanNode* node);
    bool isHierarchical() { return true; }
        
    void printPhraseTable(const WordSymbolSet & eVocab, const WordSymbolSet & fVocab, std::ostream & ptos) {
        throw std::runtime_error("HierModel::printPhraseTable not implemented");
    }

};

}

#endif
