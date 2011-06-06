
#include "pialign/parse-chart.h"

using namespace pialign;
using namespace std;

void ParseChart::initialize(int eLen, int fLen) {
#ifdef DEBUG_ON
    if(debug_)
        cerr << "Initialize: "<<eLen<<","<<fLen<<endl;
#endif
    eMultiplier_ = fLen*(fLen+1)/2+fLen+1;
    int maxSize = findChartPosition(eLen,eLen,fLen,fLen)+1;
    if((int)size() < maxSize)
        resize(maxSize);
    fill(begin(), begin()+maxSize, NEG_INFINITY);
    agendas_ = Agendas(eLen+fLen,SpanSet());
}

// trim and return an agenda
ProbSpanSet ParseChart::getTrimmedAgenda(int l, int histWidth, Prob probWidth, const LookAhead & look) { 
    SpanSet & agenda = agendas_[l-1];
    // build the beam
    ProbSpanSet ret;
    ret.reserve(agenda.size());
    for(SpanSet::const_iterator it = agenda.begin(); it != agenda.end(); it++) {
        // cerr << "Span: " << *it << " getFromChart("<<*it<<") + look.getLookAhead("<<*it<<") == "<<getFromChart(*it)<<" + "<<look.getLookAhead(*it)<<endl;
        ret.push_back(ProbSpan(getFromChart(*it)+look.getLookAhead(*it),*it));
    }
    // trim the beam if necessary
    if(probWidth != 0 || (histWidth != 0 && (int)agenda.size() > histWidth)) {
        int i, myMax = ret.size();
        if(histWidth) myMax = min(histWidth,myMax);
        sort(ret.begin(), ret.end(), moreProb);
        if(probWidth) {
            for(i = 0; i < myMax && ret[i].first>ret[0].first+probWidth; i++);
            myMax = i;
        }
        ret.resize(myMax);
    }
    return ret;
}
