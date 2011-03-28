#ifndef PROB_HIER_H__
#define PROB_HIER_H__

#include "model-base.h"
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

    void addGen(WordId jId, const Span & mySpan, Prob prob, int* tCounts, WordString & jIds) {
        phrases_.addExisting(jId);
        addAverageDerivation(jId,phrases_.getTotal(jId),prob);
    }

    void addBase(WordId jId, const Span & mySpan, Prob prob, int* tCounts, WordString & jIds) {
        phrases_.addNew(jId, -1, -1, TYPE_TERM);
        addAverageDerivation(jId,phrases_.getTotal(jId),prob);
        addType(TYPE_TERM);
    }

    void addTree(WordId jId, WordId lId, WordId rId, 
                const Span & jSpan, const Span & lSpan, const Span & rSpan, 
                int type, Prob prob, int* tCounts, WordString & jIds) {
        phrases_.addNew(jId, lId, rId, type);
        addAverageDerivation(jId,phrases_.getTotal(jId),prob);
        addType(type);
    }

    void removeSentence(WordId head, WordString & jIds, int* counts) {
        if(head == -1)
            return;
        std::vector<int> ids = phrases_.remove(head);
        // std::cerr << " phrases_.isRecursive? "<<phrases_.isRecursive()<<", removing "<<head<<",";
        for(int i = 0; i < (int)ids.size(); i++) {
            // std::cerr << " " << ids[i];
            removeType(ids[i]);
        }
        // std::cerr<<std::endl;
    }
    bool isHierarchical() { return true; }
    

    
    void printPhraseTable(const WordSymbolSet & eVocab, const WordSymbolSet & fVocab, std::ostream & ptos) {
        throw std::runtime_error("HierModel::printPhraseTable not implemented");
    }

};

}

#endif
