#ifndef BASE_COMPOUND_H__
#define BASE_COMPOUND_H__

#include "pialign/base-measure.h"
#include "pialign/dirichletdist.h"
#include <algorithm>
#include <vector>

namespace pialign {

class SpanProbElements : public SpanProbMap {

protected:

    std::vector<SpanProbMap*> maps_;

public:

    virtual ~SpanProbElements() {
        for(unsigned i = 0; i < maps_.size(); i++)
            if(maps_[i])
                delete maps_[i];
    }

    void addMap(SpanProbMap* myChart, Prob myProb) {
        maps_.push_back(myChart);
        for(SpanProbMap::const_iterator it = myChart->begin(); it != myChart->end(); it++) {
            SpanProbMap::iterator bit = this->find(it->first);
            if(bit != this->end()) {
                PRINT_DEBUG("meas span "<<it->first<<", base="<<bit->second<<", "<<myProb<<"+"<<it->second<<std::endl, 2);
                bit->second = addLogProbs(bit->second,myProb+it->second);
            } else {
                PRINT_DEBUG("meas span "<<it->first<<", nobase, "<<myProb<<"+"<<it->second<<std::endl, 2);
                this->insert(SpanProbMap::value_type(it->first,myProb+it->second));
            }
        }
    }

    std::vector<Prob> getElems(const Span & mySpan) const {
        PRINT_DEBUG("maps_.size() == "<<maps_.size()<<std::endl, 2);
        std::vector<Prob> probs(maps_.size());
        for(unsigned i = 0; i < maps_.size(); i++)
            probs[i] = maps_[i]->getProb(mySpan);
        return probs;
    }

};

class BaseCompound : public BaseMeasure {

protected:

    std::vector<BaseMeasure*> measures_;
    DirichletDist<int> dist_;

public:

    BaseCompound() : BaseMeasure(), dist_(1.0,0) { }
    ~BaseCompound() {
        for(unsigned i = 0; i < measures_.size(); i++)
            delete measures_[i];
        if(dist_.getTotal() > 0) {
            std::cerr << "WARNING: non-empty base measure " << dist_.getTotal() << std::endl;
        }
    }

    void addMeasure(BaseMeasure* meas) {
        measures_.push_back(meas);
        dist_ = DirichletDist<int>(1.0,measures_.size());
    }

    void trainPoisson(Prob avgLen, Prob nullProb) {
        for(unsigned i = 0; i < measures_.size(); i++)
            measures_[i]->trainPoisson(avgLen,nullProb);
    }

    void trainUnigrams(const Corpus & es, int eSize, const Corpus & fs, int fSize) {
        for(unsigned i = 0; i < measures_.size(); i++)
            measures_[i]->trainUnigrams(es,eSize,fs,fSize);
    }

    void setDebug(int debug) { 
        debug_ = debug;
        for(unsigned i = 0; i < measures_.size(); i++)
            measures_[i]->setDebug(debug);
    }
    void setMaxLen(int maxLen) { 
        maxLen_ = maxLen;
        for(unsigned i = 0; i < measures_.size(); i++)
            measures_[i]->setMaxLen(maxLen);
    }
    int getMaxLen() const { return maxLen_; }

    virtual SpanProbMap * getBaseChart(const WordString & e, const WordString & f) const;

    virtual void add(Span & span, WordId pid, Prob baseProb, const std::vector<Prob> & baseElems) {
        if((int)baseProbs_.size() <= pid) {
            baseProbs_.resize(pid+1,NEG_INFINITY);
            baseElems_.resize(pid+1);
        }
        baseProbs_[pid] = baseProb;
        baseElems_[pid] = baseElems;
        // find which part to add it to
        std::vector<Prob> myProbs = baseElems;
        if(baseElems.size() == 0) 
            std::cerr << "baseElems.size() == 0 for add word "<<pid << std::endl;
        else {
            PRINT_DEBUG("baseElems.size() == " <<baseElems.size()<<std::endl, 2);
            for(unsigned i = 0; i < myProbs.size(); i++)
                myProbs[i] += log(dist_.getProb(i));
            normalizeLogProbs(myProbs);
            int ans = discreteSample(myProbs,1.0);
            dist_.add(ans);
            PRINT_DEBUG("Added "<<ans<<", probability "<<dist_.getProb(ans)<<std::endl, 2);
        }
    }
    
    void remove(Span & span, WordId pid, Prob baseProb, const std::vector<Prob> & baseElems) {
        int ans;
        if(baseElems.size() == 0) {
            std::cerr << "baseElems.size() == 0 for remove word "<<pid<<", skipping"<<std::endl;
        } else {
            std::vector<Prob> myProbs = baseElems;
            Prob str = dist_.getStrength(); dist_.setStrength(0);
            for(unsigned i = 0; i < myProbs.size(); i++)
                myProbs[i] += log(dist_.getProb(i));
            normalizeLogProbs(myProbs);
            ans = discreteSample(myProbs,1.0);
            dist_.remove(ans);
            dist_.setStrength(str);
            PRINT_DEBUG("Removed "<<ans<<", probability "<<dist_.getProb(ans)<<std::endl, 2);
        }
    }

};

}

#endif
