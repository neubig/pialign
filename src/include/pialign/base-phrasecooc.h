#ifndef BASE_PHRASECOOC_H__
#define BASE_PHRASECOOC_H__

#include "pialign/base-unigram.h"
#include <algorithm>
#include <set>

namespace pialign {

class BasePhraseCooc : public BaseMeasure {

protected:

    typedef std::vector< std::pair<int,int> > SuffixArray;
    StringWordSet eSymbols_, fSymbols_;
    PairProbMap jProbs_;
    Prob condCutoff_;

public:

    BasePhraseCooc() : BaseMeasure(), condCutoff_(0.1) { }

    void substringMatrix(Corpus & corp, const WordSymbolSet & vocab, int boost, std::vector<WordString> & subs, std::vector< std::set<int> > & sents, Prob coocDisc);

    void trainCooc(Corpus & es, const WordSymbolSet & eVocab, Corpus & fs, const WordSymbolSet & fVocab, Prob coocDisc);

    SpanProbMap * getBaseChart(const WordString & e, const WordString & f) const;

    void setCondCutoff(Prob v) { condCutoff_ = v; }
    Prob getCondCutoff() { return condCutoff_; }

};

}

#endif
