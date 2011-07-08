
#include "pialign/base-unigram.h"

using namespace pialign;

SpanProbMap * BaseUnigram::getBaseChart(const WordString & e, const WordString & f) const {
    SpanProbMap * baseChart = new SpanProbMap();
    int T = e.length(), V = f.length();
    Prob eProb, fProb, noSym; //, yesSym;
    for(int s = 0; s <= T; s++) {
        eProb = 0;
        int actT = std::min(s+maxLen_,T);
        for(int t = s; t <= actT; t++) {
            if(t != s) eProb += unigrams_[e[t-1]];
            for(int u = 0; u <= V; u++) {
                fProb = 0;
                int actV = std::min(u+maxLen_,V);
                for(int v = (s==t?u+1:u); v <= actV; v++) {
                    if(u != v) fProb += unigrams_[f[v-1]];
                    Span mySpan(s,t,u,v);
                    noSym = fProb+eProb+poisProbs_[t-s]+poisProbs_[v-u];
                    // std::cerr << "BaseUnigram::calcBaseProb"<<mySpan<<" == "<<fProb<<"+"<<eProb<<"+"<<poisProbs_[t-s]<<"+"<<poisProbs_[v-u]<<" == "<<noSym<<" --> "<<yesSym<<std::endl;
                    baseChart->insertProb(mySpan,noSym);       // add to the base chart
                }
            }
        }
    }
    return baseChart;
}
