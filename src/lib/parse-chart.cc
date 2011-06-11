
#include "pialign/parse-chart.h"

using namespace pialign;
using namespace std;

void ParseChart::initialize(int eLen, int fLen) {
    PRINT_DEBUG("Initialize: "<<eLen<<","<<fLen<<endl);
    eMultiplier_ = fLen*(fLen+1)/2+fLen+1;
    int maxSize = findChartPosition(eLen,eLen,fLen,fLen)+1;
    if((int)size() < maxSize)
        resize(maxSize);
    fill(begin(), begin()+maxSize, NEG_INFINITY);
    agendas_ = Agendas(eLen+fLen,SpanVec());
}

// trim and return an agenda
ProbSpanVec ParseChart::getTrimmedAgenda(int l, Prob probWidth, const SpanSet & preserve, const LookAhead & look) { 
    SpanVec & agenda = agendas_[l-1];
    // build the beam
    ProbSpanVec ret;
    ret.reserve(agenda.size());
    Prob maxProb = NEG_INFINITY;
    for(SpanVec::const_iterator it = agenda.begin(); it != agenda.end(); it++) {
        // cerr << "Span: " << *it << " getFromChart("<<*it<<") + look.getLookAhead("<<*it<<") == "<<getFromChart(*it)<<" + "<<look.getLookAhead(*it)<<endl;
        Prob prob = getFromChart(*it)+look.getLookAhead(*it);
        maxProb = max(prob,maxProb);
        ret.push_back(ProbSpan(prob,*it));
    }
    // trim the beam if necessary
    unsigned next = 0;
    for(unsigned i = 0; i < ret.size(); i++) {
        if(ret[i].first>maxProb+probWidth || preserve.find(ret[i].second) != preserve.end())
            ret[next++] = ret[i];
    }
    ret.resize(next);
    return ret;
}
