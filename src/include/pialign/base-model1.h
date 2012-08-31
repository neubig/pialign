#ifndef BASE_MODEL1_H__
#define BASE_MODEL1_H__

#include "pialign/base-measure.h"
#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
#include <climits>

namespace pialign {

class BaseModelOne : public BaseMeasure {

protected:

    // the probability to use if the conditional probability doesn't
    //  exist (necessary when loading trimmed GIZA++ probabilities, etc.)
#define MIN_MODEL1_PROB 1e-10

    PairProbMap conds_;
    bool geometric_;

public:

    BaseModelOne() : BaseMeasure(), geometric_(false) { }

    // get the conditional probabilitiy
    Prob getCond(WordId e, WordId f) const {
        PairProbMap::const_iterator it = conds_.find(WordPairHash(e,f));
        return it != conds_.end() ? it->second : MIN_MODEL1_PROB;
    }

    // calculate combined values of P(e_i|f_{j,j+l-1}) and store them in
    // ret[i][l][j]
    void combineBases(const WordString & e, const WordString & f, std::vector<Prob> & probs) const;

    SpanProbMap * getBaseChart(const WordString & e, const WordString & f) const;
    
    void loadModelOne(const char* e2fFile, WordSymbolSet & eVocab, WordSymbolSet & fVocab, bool forward);
    
    void trainModelOne(const Corpus & es, const Corpus & fs, int eSize, int fSize, int sentLen);

    void setGeometric(bool geometric) { 
        geometric_ = geometric;
    }

    // calculate the model one probability of P(e|f)
    Prob phraseModelOne(WordString e, WordString f);

};

}

#endif
