
#include "pialign/base-compound.h"

using namespace pialign;

SpanProbMap * BaseCompound::getBaseChart(const WordString & e, const WordString & f) const {
    SpanProbElements * baseChart = new SpanProbElements();
    for(unsigned i = 0; i < measures_.size(); i++) {
        SpanProbMap * myChart = measures_[i]->getBaseChart(e,f);
        Prob myProb = log(dist_.getProb(i));
        baseChart->addMap(myChart,myProb);
    }
    return baseChart;
}
