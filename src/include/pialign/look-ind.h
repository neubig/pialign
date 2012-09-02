#ifndef LOOK_IND_H
#define LOOK_IND_H

#include "pialign/look-base.h"
#include "pialign/definitions.h"
#include <vector>
#include <algorithm>

namespace pialign {

class LookAheadInd : public LookAhead {

protected:

    std::vector<Prob> eFor_, eBack_, fFor_, fBack_, eBuff_, fBuff_;
    bool add_;
    
public:

    LookAheadInd() : add_(false) { }
    ~LookAheadInd() { }

    void setAdd(bool add) { add_ = add; }

    Prob getLookAhead(const Span & s) const {
       PRINT_DEBUG("min(eFor_[" << s.es << "]+eBack_["<<s.ee<<"], fFor_[" << s.fs << "]+fBack_["<<s.fe<<"]) == "<<
           " min("<<eFor_[s.es]<<"+"<<eBack_[s.ee]<<","<<fFor_[s.fs]<<"+"<<fBack_[s.fe]<<") == "<<
           std::min(eFor_[s.es]+eBack_[s.ee],fFor_[s.fs]+fBack_[s.fe]) << std::endl, 2);
        return std::min(eFor_[s.es]+eBack_[s.ee],fFor_[s.fs]+fBack_[s.fe]);
    }
    
    inline Prob getSpan(int start, int end, int len, const std::vector<Prob> & buff) {
        return buff[start*len+end];
    }

    inline void updateSpan(int start, int end, int len, Prob p, std::vector<Prob> & buff) {
        int idx = start*len+end;
        buff[idx] = std::max(buff[idx],p);
    }

    void addForProbs(const std::vector<Prob> & buff, std::vector<Prob> & forward, int len, int maxLen) {
        std::fill(forward.begin(),forward.end(),NEG_INFINITY);
        forward[0] = 0;
        if(!add_) {
            for(int end = 1; end <= len; end++)
                for(int start = std::max(0,end-maxLen); start < end; start++)
                    forward[end] = std::max(forward[end],forward[start]+getSpan(start,end,len,buff));
        } else {
            for(int end = 1; end <= len; end++)
                for(int start = std::max(0,end-maxLen); start < end; start++)
                    forward[end] = addLogProbs(forward[end],forward[start]+getSpan(start,end,len,buff));
        }
            
        // --- print ---
        if(GlobalVars::globalDebug >= 2) {
            PRINT_DEBUG("Forward:", 2);
            for(int i = 0; i <= len; i++)
                PRINT_DEBUG(" " << forward[i], 2);
            PRINT_DEBUG(std::endl, 2);
        }
    }
    void addBackProbs(const std::vector<Prob> & buff, std::vector<Prob> & backward, int len, int maxLen) {
        std::fill(backward.begin(),backward.end(),NEG_INFINITY);
        backward[len] = 0;
        if(!add_) {
            for(int start = len-1; start >= 0; start--)
                for(int end = start+maxLen; end > start; end--)
                    backward[start] = std::max(backward[start],backward[end]+getSpan(start,end,len,buff));
        } else {
            for(int start = len-1; start >= 0; start--)
                for(int end = start+maxLen; end > start; end--)
                    backward[start] = addLogProbs(backward[start],backward[end]+getSpan(start,end,len,buff));
        }

        // --- print ---
        if(GlobalVars::globalDebug >= 2) {
            PRINT_DEBUG("Backward:", 2);
            for(int i = 0; i <= len; i++)
                PRINT_DEBUG(" " << backward[i], 2);
            PRINT_DEBUG(std::endl, 2);
        }
    }
    
    void preCalculate(const WordString & e, const WordString & f, const SpanProbMap & base, const SpanProbMap & gen, const ParseChart & chart) {    
        // clear the buffers
        int eMax = 0, fMax = 0;
        unsigned eSize = (e.length()+1);
        if(eFor_.size() <= eSize) { eBuff_.resize(eSize*eSize); eFor_.resize(eSize); eBack_.resize(eSize); }
        std::fill( eBuff_.begin(), eBuff_.begin()+(eSize*eSize), NEG_INFINITY );
        unsigned fSize = (f.length()+1);
        if(fFor_.size() <= fSize) { fBuff_.resize(fSize*fSize); fFor_.resize(fSize); fBack_.resize(fSize); }
        std::fill( fBuff_.begin(), fBuff_.begin()+(fSize*fSize), NEG_INFINITY );
        // make the buffers
        for(SpanProbMap::const_iterator it = base.begin(); it != base.end(); it++) {
            PRINT_DEBUG("adding base "<<it->first<<", "<<it->second<<std::endl, 2);
            Prob prob = chart.getFromChart(it->first);
            eMax = std::max(it->first.ee-it->first.es,eMax);
            fMax = std::max(it->first.fe-it->first.fs,fMax);
            updateSpan(it->first.es,it->first.ee,e.length(),prob,eBuff_);
            updateSpan(it->first.fs,it->first.fe,f.length(),prob,fBuff_);
        }
        for(SpanProbMap::const_iterator it = gen.begin(); it != gen.end(); it++) {
            PRINT_DEBUG("adding gen "<<it->first<<", "<<it->second<<std::endl, 2);
            Prob prob = chart.getFromChart(it->first);
            eMax = std::max(it->first.ee-it->first.es,eMax);
            fMax = std::max(it->first.fe-it->first.fs,fMax);
            updateSpan(it->first.es,it->first.ee,e.length(),prob,eBuff_);
            updateSpan(it->first.fs,it->first.fe,f.length(),prob,fBuff_);
        }
        // count the probabilities
        addForProbs(eBuff_,eFor_,e.length(),eMax);
        addForProbs(fBuff_,fFor_,f.length(),fMax);
        addBackProbs(eBuff_,eBack_,e.length(),eMax);
        addBackProbs(fBuff_,fBack_,f.length(),fMax);
    }
    

};

}

#endif
