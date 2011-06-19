
#include "pialign/parse-chart.h"

using namespace pialign;
using namespace std;

inline void setQueue(vector< vector< pair<int,int> > > & queue, int size) {
    PRINT_DEBUG("Making queue of size "<<size<<endl);
    if((int)queue.size() <= size) queue.resize(size);
    fill(queue.begin(), queue.begin()+size, vector< pair<int,int> >());

}

void ParseChart::initialize(int eLen, int fLen) {
    PRINT_DEBUG("Initialize: "<<eLen<<","<<fLen<<endl);
    fLen_ = fLen; eLen_ = eLen;
    eMultiplier_ = fLen*(fLen+1)/2+fLen+1;
    if(useQueue_) {
        int queueSize = (fLen+1)*(eLen+1);
        setQueue(topLefts_,queueSize); setQueue(botLefts_,queueSize);
        setQueue(topRights_,queueSize); setQueue(botRights_,queueSize);
    }
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
        Prob prob = getFromChart(*it)+look.getLookAhead(*it);
        maxProb = max(prob,maxProb);
        PRINT_DEBUG("Span: " << *it << " getFromChart("<<*it<<") + look.getLookAhead("<<*it<<") == "<<getFromChart(*it)<<" + "<<look.getLookAhead(*it)<<" == "<<prob<< " (max="<<maxProb<<")"<<endl);
        ret.push_back(ProbSpan(prob,*it));
    }
    // trim the beam if necessary
    unsigned next = 0;
    PRINT_DEBUG("Max is "<<maxProb<<", saving:");
    for(unsigned i = 0; i < ret.size(); i++) {
        if(ret[i].first>maxProb+probWidth || preserve.find(ret[i].second) != preserve.end()) {
            ret[next++] = ret[i];
            PRINT_DEBUG(" "<<ret[i].first<<"/"<<ret[i].second);
        }
    }
    PRINT_DEBUG(endl);
    ret.resize(next);
    return ret;
}
