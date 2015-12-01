#ifndef PARSE_CHART_H__
#define PARSE_CHART_H__

namespace pialign {
class ParseChart;
}

#include <vector>
#include "pialign/definitions.h"
#include "pialign/look-base.h"
#include "gng/samp-gen.h"

namespace pialign {

class ParseChart : public std::vector<Prob> {

private:

    // a sorting for probability pairs
    static bool moreProb(const ProbSpan & a, const ProbSpan & b) {
        return a.first > b.first;
    }

protected:

    int eMultiplier_;
    Agendas agendas_; // the current set of bucketed agendas
    int debug_;
    bool useQueue_;
    bool viterbi_;
    int eLen_, fLen_;
    std::vector< std::vector< std::pair<int,int> > > topLefts_, botLefts_, topRights_, botRights_;

public:

    ParseChart() : std::vector<Prob>(), debug_(0), useQueue_(false), viterbi_(false) { }

    // chart access
    inline int findChartPosition(const Span & s) const {
        return findChartPosition(s.es,s.ee,s.fs,s.fe);
    }
    inline int findChartPosition(int s, int t, int u, int v) const {
        return (t*(t+1)/2+s)*eMultiplier_+v*(v+1)/2+u;
    }

    void setUseQueue(bool useQueue) { useQueue_ = useQueue; }
    void setViterbi(bool viterbi) { viterbi_ = viterbi; }

    inline std::vector< std::pair<int,int> > & getQueue(std::vector< std::vector< std::pair<int,int> > > & myQueue, int ide, int idf) {
        PRINT_DEBUG("getQueue("<<ide<<", "<<idf<<") == "<<ide*(fLen_+1)+idf<<std::endl, 2);
        return myQueue[ide*(fLen_+1)+idf];
    }
    inline std::vector< std::pair<int,int> > & getTopLefts(int ide, int idf) { return getQueue(topLefts_,ide,idf); }
    inline std::vector< std::pair<int,int> > & getBotLefts(int ide, int idf) { return getQueue(botLefts_,ide,idf); }
    inline std::vector< std::pair<int,int> > & getTopRights(int ide, int idf) { return getQueue(topRights_,ide,idf); }
    inline std::vector< std::pair<int,int> > & getBotRights(int ide, int idf) { return getQueue(botRights_,ide,idf); }
    inline void addToQueue(std::vector< std::vector< std::pair<int,int> > >& myQueue, int ide, int idf, int vale, int valf) {
        getQueue(myQueue,ide,idf).push_back(std::pair<int,int>(vale,valf));
    }

    inline Prob addToChart(int s, int t, int u, int v, Prob p, int idx, int len) {
#ifdef DEBUG_ON
        if(len == 0)
            throw std::runtime_error("ParseChart attempted to add an empty span");
        if(p != p)
            throw std::runtime_error("ParseChart attempted to add nan value");
#endif
        Prob ret;
        if((*this)[idx] <= NEG_INFINITY) {
            Span mySpan(s,t,u,v);
            agendas_[len-1].push_back(mySpan);
            if(useQueue_) {
                PRINT_DEBUG("addToQueue(topLefts_,"<<t<<","<<v<<","<<s<<","<<u<<")"<<std::endl, 2);
                PRINT_DEBUG("addToQueue(botLefts_,"<<t<<","<<u<<","<<s<<","<<v<<")"<<std::endl, 2);
                PRINT_DEBUG("addToQueue(topRights_,"<<s<<","<<v<<","<<t<<","<<u<<")"<<std::endl, 2);
                PRINT_DEBUG("addToQueue(botRights_,"<<s<<","<<u<<","<<t<<","<<v<<")"<<std::endl, 2);
                addToQueue(topLefts_,t,v,s,u);
                addToQueue(botLefts_,t,u,s,v);
                addToQueue(topRights_,s,v,t,u);
                addToQueue(botRights_,s,u,t,v);
            }
            ret = p;
        }
        else if (viterbi_)
            ret = std::max((*this)[idx], p);
        else
            ret = addLogProbs((*this)[idx],p);
        (*this)[idx] = ret;
#ifdef DEBUG_ON
        if( p <= NEG_INFINITY || p > 0 )
            throw std::runtime_error("ParseChart Attempted to add illegal probability");
        // if(debug_)
        //     std::cerr << "addToChart(Span("<<s<<","<<t<<","<<u<<","<<v<<"), "<<p<<") == " << (*this)[idx] <<std::endl;
#endif
        return ret;
    }

    inline Prob addToChart(int s, int t, int u, int v, Prob p) {
        return addToChart(s,t,u,v,p,findChartPosition(s,t,u,v),v - u + t - s);
    }
    inline Prob addToChart(const Span & s, Prob p) {
        return addToChart(s.es,s.ee,s.fs,s.fe,p);
    }

    inline Prob getFromChart(int s, int t, int u, int v) const {
        return (*this)[findChartPosition(s,t,u,v)];
    }
    inline Prob getFromChart(const Span & s) const {
        return (*this)[findChartPosition(s)];
    }
    inline void removeFromChart(int s, int t, int u, int v) {
        removeFromChart(Span(s,t,u,v));
    }
    inline void removeFromChart(const Span & s) {
        (*this)[findChartPosition(s)] = NEG_INFINITY;
    }

    // initialize the parse chart to the appropriate size
    void initialize(int eLen, int fLen);

    // trim and return an agenda
    ProbSpanVec getTrimmedAgenda(int l, Prob probWidth, const SpanSet & preserve, const LookAhead & look);

    void setDebug(int debug) { debug_ = debug; }

};

}

#endif
