#ifndef BASE_MEASURE_H__
#define BASE_MEASURE_H__

#include <vector>
#include "pialign/definitions.h"
#include "pialign/parse-chart.h"
#include "pialign/model-base.h"

namespace pialign {

class BaseMeasure {

protected:

    Prob avgLen_;
    int maxLen_;

    std::vector<Prob> poisProbs_;
    std::vector<Prob> unigrams_;

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
    virtual void addBases(const WordString & e, const WordString & f, const ProbModel & mod, ParseChart & chart, SpanProbMap & baseChart) const = 0;

    void trainPoisson(Prob avgLen, Prob nullProb) {

        // initialize poisson probabilities
        Prob poisSum = 0.0;
        poisProbs_.resize(maxLen_+1,0);
        poisProbs_[0] = exp(-1*avgLen);
        for(int i = 1; i <= maxLen_; i++) 
            poisProbs_[i] = poisProbs_[i-1] * avgLen / i;
        for(int i = 1; i <= maxLen_; i++)
            poisSum += poisProbs_[i];
        for(int i = 1; i <= maxLen_; i++)
            poisProbs_[i] = log(poisProbs_[i]/poisSum*(1-nullProb));
        poisProbs_[0] = log(nullProb); 

    }

    void trainUnigrams(const Corpus & es, int eSize, const Corpus & fs, int fSize) {
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
    void setDebug(int debug) { debug_ = debug; }
    void setMaxLen(int maxLen) { maxLen_ = maxLen; }
    int getMaxLen() const { return maxLen_; }

};

}

#endif
