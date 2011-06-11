
#include "pialign/base-unigram.h"

using namespace pialign;

void BaseUnigram::addBases(const WordString & e, const WordString & f, const ProbModel & mod, ParseChart & chart, SpanProbMap & baseChart) const {
    int T = e.length(), V = f.length();
    Prob eProb, fProb, noSym, yesSym;
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
                    noSym = uniPen_+fProb+eProb+poisProbs_[t-s]+poisProbs_[v-u];
                    yesSym = mod.calcBaseProb(mySpan, noSym);
                    // std::cerr << "BaseUnigram::calcBaseProb"<<mySpan<<" == "<<uniPen_<<"+"<<fProb<<"+"<<eProb<<"+"<<poisProbs_[t-s]<<"+"<<poisProbs_[v-u]<<" == "<<noSym<<" --> "<<yesSym<<std::endl;
                    chart.addToChart(mySpan,yesSym); // add to the overall chart
                    baseChart.insertProb(mySpan,noSym);       // add to the base chart
                }
            }
        }
    }
}
