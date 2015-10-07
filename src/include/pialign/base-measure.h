#ifndef BASE_MEASURE_H__
#define BASE_MEASURE_H__

#include <vector>
#include "pialign/definitions.h"
#include "pialign/parse-chart.h"

namespace pialign {

class BaseMeasure {

protected:

    Prob avgLen_e_;
    Prob avgLen_f_;
    int maxLen_e_;
    int maxLen_f_;

    std::vector<Prob> poisProbs_e_;
    std::vector<Prob> poisProbs_f_;
    std::vector<Prob> unigrams_;
    std::vector<Prob> baseProbs_;
    std::vector< std::vector<Prob> > baseElems_;

    // SpanProbMap baseChart_;

    int debug_;

private:

    int addSums(const Corpus & corp, std::vector<int> & sums) {
        int ret = 0;
        for(int i = 0; i < (int)corp.size(); i++) {
            ret += corp[i].length();
            for(int j = 0; j < (int)corp[i].length(); j++)
                sums[corp[i][j]]++;
        }
        return ret;
    }

public:

    BaseMeasure() : debug_(0) { }

    virtual ~BaseMeasure() { }
    
    // add base probabilities for the unigram model
    //  e,f strings, logPenalty is the additional log probability added
    //  by fallbacks, etc in the model, chart is the overall chart
    virtual SpanProbMap * getBaseChart(const WordString & e, const WordString & f) const = 0;

    virtual void trainPoisson(Prob avgLenE, Prob nullProbE, Prob avgLenF, Prob nullProbF) {

        // initialize poisson probabilities (E side)
        Prob poisSum = 0.0;
        poisProbs_e_.resize(maxLen_e_+1,0);
        poisProbs_e_[0] = exp(-1*avgLenE);
        for(int i = 1; i <= maxLen_e_; i++) 
            poisProbs_e_[i] = poisProbs_e_[i-1] * avgLenE / i;
        for(int i = 1; i <= maxLen_e_; i++)
            poisSum += poisProbs_e_[i];
        for(int i = 1; i <= maxLen_e_; i++)
            poisProbs_e_[i] = log(poisProbs_e_[i]/poisSum*(1-nullProbE));
        poisProbs_e_[0] = log(nullProbE); 

        // initialize poisson probabilities (F side)
        poisSum = 0.0;
        poisProbs_f_.resize(maxLen_f_+1,0);
        poisProbs_f_[0] = exp(-1*avgLenF);
        for(int i = 1; i <= maxLen_f_; i++) 
            poisProbs_f_[i] = poisProbs_f_[i-1] * avgLenF / i;
        for(int i = 1; i <= maxLen_f_; i++)
            poisSum += poisProbs_f_[i];
        for(int i = 1; i <= maxLen_f_; i++)
            poisProbs_f_[i] = log(poisProbs_f_[i]/poisSum*(1-nullProbF));
        poisProbs_f_[0] = log(nullProbF); 
    }

    virtual void trainUnigrams(const Corpus & es, int eSize, const Corpus & fs, int fSize) {
        // make the unigram probabilities
        std::vector<int> sums(eSize+fSize,0);
        unigrams_.resize(sums.size(),0);
        double eTot = addSums(es,sums);
        double fTot = addSums(fs,sums);
        for(int i = 0; i < eSize; i++) unigrams_[i] = log(sums[i]/eTot);
        for(int i = eSize; i < eSize+fSize; i++) unigrams_[i] = log(sums[i]/fTot);
    }

    void initialize() { }
//     void addToChart(const Span & span, Prob myProb) {
//         baseChart_.insert(std::pair<Span,Prob>(span,myProb));
//     }
//    Prob getFromChart(const Span & mySpan) {
//        SpanProbMap::const_iterator it = baseChart_.find(mySpan);
//        return it == baseChart_.end() ? NEG_INFINITY : it->second;
//
//    }
    virtual void setDebug(int debug) { debug_ = debug; }
    virtual void setMaxLen(int maxLenE, int maxLenF) { maxLen_e_ = maxLenE; maxLen_f_ = maxLenF; }
    virtual int getMaxLenE() const { return maxLen_e_; }
    virtual int getMaxLenF() const { return maxLen_f_; }

    virtual void add(Span & span, WordId pid, Prob baseProb, const std::vector<Prob> & baseElems) {
        if((int)baseProbs_.size() <= pid) {
            baseProbs_.resize(pid+1,NEG_INFINITY);
            baseElems_.resize(pid+1);
        }
        baseProbs_[pid] = baseProb;
        baseElems_[pid] = baseElems;
    }

    virtual void remove(Span & span, WordId pid, Prob baseProb, const std::vector<Prob> & baseElems) { }

    virtual Prob getBase(WordId pid) {
        if((int)baseProbs_.size() <= pid)
            THROW_ERROR("Phrase id "<<pid<<" is larger than size "<<baseProbs_.size());
        return baseProbs_[pid];
    }
    
    virtual std::vector<Prob> getElems(WordId pid) {
        if((int)baseElems_.size() <= pid)
            THROW_ERROR("Phrase id "<<pid<<" is larger than size "<<baseElems_.size());
        return baseElems_[pid];
    }

};

}

#endif
