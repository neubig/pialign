#ifndef LOOK_BASE_H
#define LOOK_BASE_H

namespace pialign {
class LookAhead;
}

#include "pialign/definitions.h"
#include "pialign/parse-chart.h"

namespace pialign {

class LookAhead {
    
public:

    LookAhead() { }
    virtual ~LookAhead() { }

    virtual Prob getLookAhead(const Span & s) const = 0;

    virtual void preCalculate(const WordString & e, const WordString & f, const SpanProbMap & base, const SpanProbMap & gen, const ParseChart & chart) = 0;    

};

}

#endif
