#ifndef LOOK_BASE_H
#define LOOK_BASE_H

#include "pialign/definitions.h"

namespace pialign {

class LookAhead {
    
public:

    LookAhead() { }
    virtual ~LookAhead() { }

    virtual Prob getLookAhead(const Span & s) const = 0;

    virtual void preCalculate(const WordString & e, const WordString & f, const SpanProbMap & baseChart, const SpanProbMap & genChart) = 0;    

};

}

#endif
