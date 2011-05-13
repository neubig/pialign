#ifndef DIRICHLETDIST_H__
#define DIRICHLETDIST_H__

#include "pialign/definitions.h"
#include "gng/samp-gen.h"
#include <vector>
#include <stdexcept>
#include <cmath>
#include <iostream>

namespace pialign {

template < class T >
class DirichletDist {

private:

    T total_;
    std::vector< T > counts_;
    std::vector< Prob > bases_;

    Prob stren_;

public:

    DirichletDist(Prob stren, T size) : 
        total_(0), counts_(size,0), bases_(size, 1.0/size), stren_(stren) { }

    Prob getProb(T id) const {
        return (counts_[id]+bases_[id]*stren_)/(total_+stren_);
    }

    Prob getFallbackProb() const {
        return stren_/(total_+stren_);
    }

    Prob getLogProb(WordId id) const {
        return log(getProb(id));
    }

    std::vector<Prob> getAllProbs() const {
        std::vector<Prob> ret();
        for(T i = 0; i < ret.size(); i++)
            ret[i] = getProb(i);
        return ret;
    }

    std::vector<Prob> getAllLogProbs(const std::vector<Prob> & bases) const {
        std::vector<Prob> ret = getAllProbs(bases);
        for(T i = 0; i < ret.size(); i++)
            ret[i] = log(ret[i]);
        return ret;
    }

    void add(T id) {
        counts_[id]++;
        total_++;
    }
     
    void remove(T id) {
        counts_[id]--;
        total_--;
#ifdef DEBUG_ON
        if(counts_[id] < 0)
            throw std::runtime_error("Count less than zero for dirichlet dist");
#endif
    }

    Prob getStrength() const { return stren_; }
    void setStrength(Prob stren) { stren_ = stren; }
    T getTotal() const { return total_; }
    int getSize() const { return counts_.size(); }
    const std::vector<T> & getCounts() const { return counts_; }
    void setBases(const std::vector< Prob > & bases) { bases_ = bases; }

};

}

#endif
