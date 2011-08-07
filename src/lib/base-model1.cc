#include "pialign/base-model1.h"

using namespace pialign;
using namespace std;

void BaseModelOne::combineBases(const WordString & e, const WordString & f, std::vector<Prob> & probs) const {
    int I = e.length(), J = f.length()+1, myMax = I*J*(maxLen_+1);
    if((int)probs.size() < myMax) probs.resize(myMax);
    fill(probs.begin(),probs.begin()+myMax,0);
    for(int i = 0; i < I; i++) {
        Prob* base = &probs[i*J*(maxLen_+1)+J];
        int myE = e[i];
        Prob nullProb = getCond(myE,0);
        for(int j = 0; j < J-1; j++) {
            // add zero and one values
            base[j-J] = nullProb;
            base[j] = getCond(myE,f[j]);
            PRINT_DEBUG("combineBases("<<i<<","<<j<<") == "<<base[j]<<std::endl)
        }
        base[-1] = nullProb;
        for(int l = 1; l < maxLen_; l++)
            for(int j = 0; j < J-l-1; j++) {
                PRINT_DEBUG("combineBases("<<i<<","<<j<<","<<l<<") == "<<base[j]<<std::endl);
                base[j+l*J] = (base[j+(l-1)*J]*l+base[j+l])/(l+1);
            }
    }
    for(int i = 0; i < myMax; i++) {
        probs[i] = probs[i] ? log(probs[i]) : NEG_INFINITY;
    }
}

SpanProbMap * BaseModelOne::getBaseChart(const WordString & e, const WordString & f) const {
    SpanProbMap * baseChart = new SpanProbMap();

    int T = e.length(), V = f.length();
    Prob l2 = log(2);
    Prob eProb, fProb, noSym, fm1;
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
                        fm1 += fComb[(v-1)*(T+1)*(maxLen_+1)+(t-s)*(T+1)+s];
                    }
                    Prob em1 = 0;
                    if(s != t) 
                        em1 = (em1s[v*(maxLen_+1)+v-u] += eComb[(t-1)*(V+1)*(maxLen_+1)+(v-u)*(V+1)+u]);
                    Span mySpan(s,t,u,v);
                    // add model one probabilities 
                    if(geometric_) 
                        noSym = (em1+fProb+fm1+eProb)/2+poisProbs_[t-s]+poisProbs_[v-u];
                    else 
                        noSym = addLogProbs(em1+fProb,fm1+eProb)-l2+poisProbs_[t-s]+poisProbs_[v-u];
                    // PRINT_DEBUG("calcBaseProb @ "<<mySpan<<", addLogProbs("<<em1<<"+"<<fProb<<","<<fm1<<"+"<<eProb<<")+"<<poisProbs_[t-s]<<"+"<<poisProbs_[v-u]<<") == "<<noSym<<" --> "<<yesSym<<std::endl);
                    baseChart->insertProb(mySpan,noSym); // add to the base chart
                }
            }
        }
    }
    return baseChart;
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
    int maxIters = 100;
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
                }
            }
        }
        // M step
        for(PairProbMap::iterator it = count.begin(); it != count.end(); it++)
            conds_[it->first] = it->second/total[it->first.second];
        std::cerr << " Iteration " << ++iter << ": likelihood "<<lik<<std::endl;
    } while((lastLik == 0.0 || (lastLik-lik) < lik * likCut) && --maxIters > 0);
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
