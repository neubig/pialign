#ifndef PYDIST_H__
#define PYDIST_H__

#include "pialign/definitions.h"
#include "gng/samp-gen.h"
#include <vector>
#include <stdexcept>
#include <cmath>
#include <iostream>

#define PRIOR_SA 2.0
#define PRIOR_SB 1.0
#define PRIOR_DA 2.0
#define PRIOR_DB 2.0

namespace pialign {

template < class T >
class PyTable
{
public:
    PyTable() : count(0), left(-1), right(-1), type(-1) { }
    PyTable(int c, T l, T r, T t) : count(c), left(l), right(r), type(t) { }
    int count;    
    T left, right, type;
};

template < class T >
class PyTableSet : public std::vector< PyTable<T> > {

public:
    int total;
    PyTableSet() : std::vector< PyTable<T> >(), total(0) { }

};

template < class T >
class PyDenseIndex {
protected:
    typedef PyTableSet<T> TSet;
    std::vector< TSet > idx_;
public:
    typedef typename std::vector< TSet >::iterator iterator;
    iterator begin() { return idx_.begin(); }
    iterator end() { return idx_.end(); }
    TSet & iterTableSet(iterator it) { return *it; }
    
    const PyTableSet<T> * getTableSet(T id) const { 
        return (int)idx_.size() > id ? &(idx_[id]) : 0; 
    }
    int getTotal(T id) const {
        return (int)idx_.size() > id ? (idx_[id]).total : 0; 
    }
    PyTableSet<T> & addTableSet(T id) {
        // std::cerr << "addTableSet("<<id<<"), size="<<idx_.size()<<std::endl;
        if((int)idx_.size() <= id) idx_.resize(id+1);
        return idx_[id];
    }
};
template < class T >
class PySparseIndex {
protected:
    typedef PyTableSet<T> TSet;
    std::tr1::unordered_map< int, TSet > idx_;
public:
    typedef typename std::tr1::unordered_map< int, TSet >::const_iterator const_iterator;
    typedef typename std::tr1::unordered_map< int, TSet >::iterator iterator;
    iterator begin() { return idx_.begin(); }
    iterator end() { return idx_.end(); }
    TSet & iterTableSet(iterator it) { return it->second; }

    const TSet * getTableSet(T id) const {
        const_iterator it = idx_.find(id);
        return it == idx_.end() ? 0 : & it->second;
    }
    int getTotal(T id) const {
        const_iterator it = idx_.find(id);
        return it == idx_.end() ? 0 : it->second.total;
    }
    TSet & addTableSet(T id) {
        iterator it = idx_.find(id);
        if(it == idx_.end()) 
            it = idx_.insert( std::pair< T, TSet >(id,TSet()) ).first;
        return it->second;
    }
};

template < class T, class Index >
class PyDist {

private:

    typedef PyTableSet<T> TSet;
    typedef typename TSet::iterator TSetIter;
    typedef typename TSet::const_iterator TSetCIter;

    T total_, tables_;
    Index counts_;
    Prob spAlpha_, spBeta_, dpAlpha_, dpBeta_;

    // whether to delete tables recursively
    bool removedTable_;

    Prob stren_, disc_;

    // the last table that was removed
    PyTable<T> lastTable_;

public:

    PyDist(Prob stren, Prob disc) : total_(0), tables_(0), counts_(), 
        spAlpha_(PRIOR_SA), spBeta_(PRIOR_SB), dpAlpha_(PRIOR_DA), dpBeta_(PRIOR_DB),
        removedTable_(false), stren_(stren), disc_(disc) { }

    Prob getProb(T id, Prob base) const {
        const TSet * tab = counts_.getTableSet(id);
        Prob myCount = tab ? tab->total-disc_*tab->size() : 0;
        return ( myCount+base*(stren_+tables_*disc_) )/(total_+stren_);
    }
    int getTotal(T id) const {
        return counts_.getTotal(id);
    }

    Prob getFallbackProb() const {
        // std::cerr << "("<<stren_<<"+"<<tables_<<"*"<<disc_<<")/("<<total_<<"+"<<stren_<<") == "<<(stren_+tables_*disc_)/(total_+stren_)<<std::endl;
        return (stren_+tables_*disc_)/(total_+stren_);
    }

    Prob getLogProb(WordId id, Prob base) const {
        return log(getProb(id,exp(base)));
    }

    std::vector<Prob> getAllProbs(const std::vector<Prob> & bases) const {
#ifdef DEBUG_ON
        if(bases.size() != counts_.size())
            throw std::runtime_error("Mismatched bases and probs");
#endif
        std::vector<Prob> ret(bases.size());
        for(T i = 0; i < ret.size(); i++)
            ret[i] = getProb(i,bases[i]);
        return ret;
    }

    std::vector<Prob> getAllLogProbs(const std::vector<Prob> & bases) const {
        std::vector<Prob> ret = getAllProbs(bases);
        for(T i = 0; i < ret.size(); i++)
            ret[i] = log(ret[i]);
        return ret;
    }

    // addExisting works when tables have been deleted and the total table
    // count is zero, but there is still one cached
    void addExisting(T id) {
        TSet & set = counts_.addTableSet(id);
#ifdef DEBUG_ON
        if(set.size() == 0)
            throw std::runtime_error("PyDist::addExisting with no tables");
#endif
        // std::cerr << " pydist::addExisting("<<id<<")"<<std::endl;
        int mySize = set.size();
        TSetIter it = set.begin();
        if(mySize > 1) {
            Prob left = rand()*(set.total-mySize*disc_)/RAND_MAX;
            while(it != set.end() && (left -= it->count-disc_) > 0)
                it++;
            if(it == set.end()) it--;
        }
        it->count++;
        set.total++;
        total_++;
    }

    void addNew(T id, T left = -1, T right = -1, T type = -1) {
        // std::cerr << " pydist::addNew("<<id<<","<<left<<","<<right<<","<<type<<")"<<std::endl;
        if(type != TYPE_TERM && type != -1 && (left == -1 || right == -1))
            throw std::runtime_error("No left and right for recursive definition");
        if(id == left || id == right)
            throw std::runtime_error("Recursive definition in addNew");
        TSet & set = counts_.addTableSet(id);
        // std::cerr << " set.size() == "<<set.size()<<std::endl;
        set.push_back(PyTable<T>(1,left,right,type));
        set.total++;
        tables_++;
        total_++;
    }

    void add(T id, Prob base, T left = -1, T right = -1, T type = -1) {
        Prob genProb = getProb(id,0), baseProb = base*getFallbackProb();
        if(bernoulliSample(genProb/(genProb+baseProb)))
            addExisting(id);
        else
            addNew(id,left,right,type);
    }
    
    // remove one customer, return the probability of
    //  the cache if removed from the cache
    //  the fallback if removed from the base
    Prob remove(T id) {
        if(id < 0) return 0;
#ifdef DEBUG_ON
        if(!counts_.getTableSet(id))
            throw std::runtime_error("Overflow in PyDist::remove"); 
#endif
        removedTable_ = false;
        Prob prob = 0;
        TSet & set = counts_.addTableSet(id);
        int mySize = set.size();
        TSetIter it = set.begin(); 
        if(mySize > 1) {
            int left = rand() % set.total;
            while((left -= it->count) >= 0)
                it++;
        }
        it->count--;
        set.total--;
        total_--;
        if(it->count == 0) {
            removedTable_ = true;
            lastTable_ = *it;
            PyTable<T> save = lastTable_;
            tables_--;
            set.erase(it);
            prob = log(getFallbackProb());
        } else 
            prob = log(getProb(id,0));
        return prob;
    }

    // sample the parameters
    void sampleParameters() {
        // std::cerr << "sampling parameters" << std::endl;
        Prob x = total_ > 1 ? betaSample(stren_+1,total_-1) : 1;
        Prob y = 0, z = 0;
        for(int i = 1; i < tables_; i++)
            y += bernoulliSample(stren_/(stren_+disc_*i));
        for(typename Index::iterator it2 = counts_.begin(); it2 != counts_.end(); it2++) {
            const PyTableSet<T> & set = counts_.iterTableSet(it2);
            for(typename PyTableSet<T>::const_iterator it = set.begin();
                    it != set.end(); it++) { 
                for(int k = 1; k < it->count; k++) {
                    z += (1-bernoulliSample((k-1)/(k-disc_)));
                }
            }
        }
        // std::cerr << "disc_ = betaSample("<<dpAlpha_<<"+"<<(tables_-1)<<"-"<<y<<","<<dpBeta_<<"+"<<z<<")"<<std::endl;
        disc_ = betaSample(dpAlpha_+(tables_-1)-y,dpBeta_+z);
        // std::cerr << "stren_ = gammaSample("<<spAlpha_<<"+"<<y<<","<<spBeta_<<"-log("<<x<<"))"<<std::endl;
        stren_ = gammaSample(spAlpha_+y,spBeta_-log(x));
        // std::cerr << "Done!" << std::endl;
    }

    Prob getStrength() const { return stren_; }
    void setStrength(Prob stren) { stren_ = stren; }
    Prob getDiscount() const { return disc_; }
    void setDiscount(Prob disc) { disc_ = disc; }
    int getTableCount() const { return tables_; }

    void setDiscountAlpha(Prob v) { dpAlpha_ = v; }
    void setDiscountBeta(Prob v) { dpBeta_ = v; }
    void setStrengthAlpha(Prob v) { spAlpha_ = v; }
    void setStrengthBeta(Prob v) { spBeta_ = v; }

    const PyTable<T> & getLastTable() const { return lastTable_; }
    bool isRemovedTable() { return removedTable_; }

};

}

#endif
