#ifndef BASE_UNIGRAM_H__
#define BASE_UNIGRAM_H__

#include "pialign/base-measure.h"
#include <algorithm>

namespace pialign {

class BaseUnigram : public BaseMeasure {

protected:

    Prob uniPen_;

public:

    BaseUnigram() : BaseMeasure(), uniPen_(0) { }

    virtual void addBases(const WordString & e, const WordString & f, const ProbModel & mod, ParseChart & chart, SpanProbMap & baseChart) const;

};

}

#endif
