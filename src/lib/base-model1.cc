#include "pialign/base-model1.h"

using namespace pialign;
using namespace std;

void BaseModelOne::combineBases(const WordString & e, const WordString & f, std::vector<Prob> & probs) const {
    int I = e.length(), J = f.length(), myMax = I*J*(maxLen_+1);
    if((int)probs.size() < myMax) probs.resize(myMax);
    for(int i = 0; i < I; i++) {
        Prob* base = &probs[i*J*(maxLen_+1)+J];
        int myE = e[i];
        Prob nullProb = getCond(myE,0);
        for(int j = 0; j < J; j++) {
            // add zero and one values values
            base[j-J] = nullProb;
            base[j] = getCond(myE,f[j]);
#ifdef DEBUG_ON
            if(debug_)
                std::cerr << "combineBases("<<i<<","<<j<<") == "<<base[j]<<std::endl;
#endif
        }
        for(int l = 1; l < maxLen_; l++)
            for(int j = 0; j < J-l; j++) {
#ifdef DEBUG_ON
                if(debug_)
                    std::cerr << "combineBases("<<i<<","<<j<<","<<l<<") == ("<<base[j+(l-1)*J]<<"*"<<l<<"+"<<base[j+l]<<")/"<<(l+1)<<std::endl;
#endif
                base[j+l*J] = (base[j+(l-1)*J]*l+base[j+l])/(l+1);
            }
    }
    for(int i = 0; i < myMax; i++)
        probs[i] = probs[i] ? log(probs[i]) : NEG_INFINITY;
}

void BaseModelOne::addBases(const WordString & e, const WordString & f, const ProbModel & mod, ParseChart & chart, SpanProbMap & baseChart) const {
    int T = e.length(), V = f.length();
    Prob l2 = log(2);
    Prob eProb, fProb, myProb, fm1;
    std::vector<Prob> em1s((V+1)*(maxLen_+1));
    std::vector<Prob> eComb, fComb;
    combineBases(e,f,eComb); combineBases(f,e,fComb);
    for(int s = 0; s <= T; s++) {
        eProb = 0;
        fill(em1s.begin(),em1s.end(),0.0);
        int actT = std::min(s+maxLen_,T);
        for(int t = s; t <= actT; t++) {
            if(t != s) eProb += unigrams_[e[t-1]];
            for(int u = 0; u <= V; u++) {
                fProb = 0; fm1 = 0;
                int actV = std::min(u+maxLen_,V);
                for(int v = (s==t?u+1:u); v <= actV; v++) {
                    if(u != v) {
                        fProb += unigrams_[f[v-1]];
                        fm1 += fComb[(v-1)*T*(maxLen_+1)+(t-s)*T+s];
                    }
                    Prob em1 = 0;
                    if(s != t) 
                        em1 = (em1s[v*(maxLen_+1)+v-u] += eComb[(t-1)*V*(maxLen_+1)+(v-u)*V+u]);
                    Span mySpan(s,t,u,v);
                    // add model one probabilities 
                    if(geometric_) {
                        myProb = mod.calcBaseProb(mySpan,(em1+fProb+fm1+eProb)/2+poisProbs_[t-s]+poisProbs_[v-u]);
                    }
                    else {
                        myProb = mod.calcBaseProb(mySpan,addLogProbs(em1+fProb,fm1+eProb)-l2+poisProbs_[t-s]+poisProbs_[v-u]);
                    }
#ifdef DEBUG_ON
                    if(debug_)
                        std::cerr << "calcBaseProb @ "<<mySpan<<", addLogProbs("<<em1<<"+"<<fProb<<","<<fm1<<"+"<<eProb<<")+"<<poisProbs_[t-s]<<"+"<<poisProbs_[v-u]<<") == "<<myProb<<std::endl;
#endif
                    chart.addToChart(mySpan,myProb); // add to the overall chart
                    baseChart.insertProb(mySpan,myProb);       // add to the base chart
                }
            }
        }
    }
}

void BaseModelOne::loadModelOne(const char* e2fFile, WordSymbolSet & eVocab, WordSymbolSet & fVocab, bool forward) {
    std::cerr << "Loading model one from "<<e2fFile<<std::endl;
    std::ifstream ifs(e2fFile);
    if(!ifs) {
        std::ostringstream oss;
        oss << "IO Error: Couldn't open "<<e2fFile;
        throw std::runtime_error(oss.str());
    }
    int eBonus = (forward?0:fVocab.size());
    int fBonus = (forward?eVocab.size():0);
    int lineCount = 0;
    std::string line,e,f;
    Prob c;
    while(getline(ifs, line)) {
        std::istringstream iss(line);
        iss >> f >> e >> c;
        WordId fId = (f == "NULL"?0:fVocab.getId(f,false));
        WordId eId = (e == "NULL"?0:eVocab.getId(e,false));
        if(eId == -1 || fId == -1) {
            std::cerr <<"warning, word in lexicon file that doesn't exist in corpus ("<<eId<<","<<fId<<") "<<line<<std::endl;
        } else if (fId != 0) {
            if(eId) eId+=eBonus;
            fId+=fBonus;
            conds_.insert(std::pair< std::pair<WordId, WordId>, Prob>(std::pair<WordId,WordId>(fId,eId), c));
            lineCount++;
        }
    }
    std::cerr << " "<<lineCount<<" lines"<<std::endl;
}

void BaseModelOne::trainModelOne(const Corpus & es, const Corpus & fs, int eSize, int fSize, int sentLen) {
    PairProbMap count;
    std::vector<Prob> total(fSize+eSize), tBuff(sentLen+1); 
    std::vector< std::pair<WordId,WordId> > pBuff(sentLen+1);
    // initialize the probability
    std::cerr << "Training model 1" << std::endl;
    int i,j,k,iter=0;
    Prob uniProb = 1.0/eSize;
    for(i = 0; i < (int)es.size(); i++) {
        // cerr << "i="<<i<<", el="<<es[i].length()<<", fl="<<fs[i].length()<<endl;
        // cerr << "i="<<i<<", j="<<j<<", k="<<k<<endl;
        for(j = 0; j < (int)es[i].length(); j++) {
            for(k = 0; k < (int)fs[i].length(); k++) {
                std::pair<WordId,WordId> id(es[i][j],fs[i][k]);
                conds_.insert(PairProbMap::value_type(id,uniProb));
                count.insert(PairProbMap::value_type(id,0.0));
            }
            std::pair<WordId,WordId> id(es[i][j],0);
            conds_.insert(PairProbMap::value_type(id,uniProb));
            count.insert(PairProbMap::value_type(id,0.0));
        }
    }
    // train the model
    Prob lastLik = 0.0, likCut = 0.001, sTotal, lik = 0.0, norm;
    do {
        // reset the values
        lastLik = lik;
        lik = 0.0;
        for(PairProbMap::iterator it = count.begin(); it != count.end(); it++)
            it->second = 0.0;
        fill(total.begin(),total.end(),0.0);
        // E step
        for(i = 0; i < (int)es.size(); i++) {
            if(es[i].length()*fs[i].length() == 0)
                continue;
            // cerr << "i="<<i<<", el="<<es[i].length()<<", fl="<<fs[i].length()<<endl;
            for(j = 0; j < (int)es[i].length(); j++) {
                sTotal = 0;
                // do words + null
                for(k = 0; k < (int)fs[i].length(); k++) {
                    pBuff[k] = std::pair<WordId,WordId>(es[i][j],fs[i][k]);
                    tBuff[k] = conds_.find(pBuff[k])->second;
                    sTotal += tBuff[k];
                }
                pBuff[k] = std::pair<WordId,WordId>(es[i][j],0);
                tBuff[k] = conds_.find(pBuff[k])->second;
                sTotal += tBuff[k];
                // likelihood
                lik += log(sTotal/(fs[i].length()+1));
                // do words + null
                for(k = 0; k < (int)fs[i].length()+1; k++) {
                    norm = tBuff[k]/sTotal;
                    count[pBuff[k]] += norm;
                    total[pBuff[k].second] += norm;
                    // std::cerr << "count["<<pBuff[k].first<<","<<pBuff[k].second<<"] += "<<norm<<" == "<<count[pBuff[k]]<<std::endl;
                    // std::cerr << "total["<<pBuff[k].second<<"] += "<<norm<<" == "<<total[pBuff[k].second]<<std::endl;
                }
            }
        }
        // M step
        for(PairProbMap::iterator it = count.begin(); it != count.end(); it++) {
            conds_[it->first] = it->second/total[it->first.second];
            // std::cerr << "conds_["<<it->first.first<<","<<it->first.second<<"] == "<<it->second<<"/"<<total[it->first.second]<<" == "<<conds_[it->first]<<std::endl;
        }
        std::cerr << " Iteration " << ++iter << ": likelihood "<<lik<<std::endl;
    } while(lastLik == 0.0 || (lastLik-lik) < lik * likCut);
}

// calculate the model one probability of P(e|f)
Prob BaseModelOne::phraseModelOne(WordString e, WordString f) {
    Prob ret = 1;
    int el = e.length(), fl = f.length();
    for(int i = 0; i < el; i++) {
        Prob myProb = getCond(e[i],0);
        for(int j = 0; j < fl; j++)
            myProb += getCond(e[i],f[j]);
        ret *= (myProb/(el+1));
    }
    return ret;
}
