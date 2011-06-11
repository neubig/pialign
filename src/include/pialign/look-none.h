#ifndef LOOK_NONE_H
#define LOOK_NONE_H

#include "pialign/definitions.h"
#include "pialign/look-base.h"

namespace pialign {

class LookAheadNone : public LookAhead {
    
public:

    LookAheadNone() { }
    ~LookAheadNone() { }

    Prob getLookAhead(const Span & s) const { return 0; }

    void preCalculate(const WordString & e, const WordString & f, const SpanProbMap & base, const SpanProbMap & gen, const ParseChart & chart) { }
    

};

}

#endif
