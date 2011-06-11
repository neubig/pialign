#ifndef BASE_PHRASECOOC_H__
#define BASE_PHRASECOOC_H__

#include "pialign/base-unigram.h"
#include <algorithm>
#include <set>

namespace pialign {

class BasePhraseCooc : public BaseUnigram {

protected:

    typedef std::vector< std::pair<int,int> > SuffixArray;
    StringWordSet eSymbols_, fSymbols_;
    PairProbMap jProbs_;

public:

    BasePhraseCooc() : BaseUnigram() { }

    void substringMatrix(Corpus & corp, const WordSymbolSet & vocab, int boost, std::vector<WordString> & subs, std::vector< std::set<int> > & sents, Prob coocDisc);

    void trainCooc(Corpus & es, const WordSymbolSet & eVocab, Corpus & fs, const WordSymbolSet & fVocab, Prob coocDisc, Prob coocCut);

    void addBases(const WordString & e, const WordString & f, const ProbModel & mod, ParseChart & chart, SpanProbMap & baseChart) const;

};

}

#endif
