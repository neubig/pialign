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
       // std::cerr << "min(eFor_[" << s.es << "]+eBack_["<<s.ee<<"], fFor_[" << s.fs << "]+fBack_["<<s.fe<<"]) == "<<
       //     " min("<<eFor_[s.es]<<"+"<<eBack_[s.ee]<<","<<fFor_[s.fs]<<"+"<<fBack_[s.fe]<<") == "<<
       //     std::min(eFor_[s.es]+eBack_[s.ee],fFor_[s.fs]+fBack_[s.fe]) << std::endl;
        return std::min(eFor_[s.es]+eBack_[s.ee],fFor_[s.fs]+fBack_[s.fe]);
    }
    
    inline Prob getSpan(int start, int end, int len, const std::vector<Prob> & buff) {
        return buff[start*len+end];
    }

    inline void updateSpan(int start, int end, int len, Prob p, std::vector<Prob> & buff) {
        int idx = start*len+end;
        buff[idx] = std::max(buff[idx],p);
    }

    void addForProbs(const std::vector<Prob> & buff, std::vector<Prob> & forward, int len) {
        std::fill(forward.begin(),forward.end(),NEG_INFINITY);
        forward[0] = 0;
        if(!add_) {
            for(int end = 1; end <= len; end++)
                for(int start = 0; start < end; start++)
                    forward[end] = std::max(forward[end],forward[start]+getSpan(start,end,len,buff));
        } else {
            for(int end = 1; end <= len; end++)
                for(int start = 0; start < end; start++)
                    forward[end] = addLogProbs(forward[end],forward[start]+getSpan(start,end,len,buff));
        }
            
        // // --- print ---
        // std::cerr << "Forward:";
        // for(int i = 0; i <= len; i++)
        //     std::cerr << " " << forward[i];
        // std::cerr << std::endl;
    }
    void addBackProbs(const std::vector<Prob> & buff, std::vector<Prob> & backward, int len) {
        std::fill(backward.begin(),backward.end(),NEG_INFINITY);
        backward[len] = 0;
        if(!add_) {
            for(int start = len-1; start >= 0; start--)
                for(int end = start+1; end <= len; end++)
                    backward[start] = std::max(backward[start],backward[end]+getSpan(start,end,len,buff));
        } else {
            for(int start = len-1; start >= 0; start--)
                for(int end = start+1; end <= len; end++)
                    backward[start] = addLogProbs(backward[start],backward[end]+getSpan(start,end,len,buff));
        }

        // // --- print ---
        // std::cerr << "Backward:";
        // for(int i = 0; i <= len; i++)
        //     std::cerr << " " << backward[i];
        // std::cerr << std::endl;
    }

    void preCalculate(const WordString & e, const WordString & f, const SpanProbMap & baseChart, const SpanProbMap & genChart) {
        // clear the buffers
        unsigned eSize = (e.length()+1);
        if(eFor_.size() <= eSize) { eBuff_.resize(eSize*eSize); eFor_.resize(eSize); eBack_.resize(eSize); }
        std::fill( eBuff_.begin(), eBuff_.begin()+(eSize*eSize), NEG_INFINITY );
        unsigned fSize = (f.length()+1)*(f.length()+1);
        if(fFor_.size() <= fSize) { fBuff_.resize(fSize*fSize); fFor_.resize(fSize); fBack_.resize(fSize); }
        std::fill( fBuff_.begin(), fBuff_.begin()+(fSize*fSize), NEG_INFINITY );
        // make the buffers
        for(SpanProbMap::const_iterator it = baseChart.begin(); it != baseChart.end(); it++) {
            // std::cerr << "adding base "<<it->first<<", "<<it->second<<std::endl;
            updateSpan(it->first.es,it->first.ee,e.length(),it->second,eBuff_);
            updateSpan(it->first.fs,it->first.fe,f.length(),it->second,fBuff_);
        }
        for(SpanProbMap::const_iterator it = genChart.begin(); it != genChart.end(); it++) {
            // std::cerr << "adding gen "<<it->first<<", "<<it->second<<std::endl;
            updateSpan(it->first.es,it->first.ee,e.length(),it->second,eBuff_);
            updateSpan(it->first.fs,it->first.fe,f.length(),it->second,fBuff_);
        }
        // count the probabilities
        addForProbs(eBuff_,eFor_,e.length());
        addForProbs(fBuff_,fFor_,f.length());
        addBackProbs(eBuff_,eBack_,e.length());
        addBackProbs(fBuff_,fBack_,f.length());
    }
    

};

}

#endif
