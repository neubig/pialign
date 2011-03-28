#ifndef PROB_LENGTH_H__
#define PROB_LENGTH_H__

#include "model-base.h"
#include <stdexcept>

namespace pialign {

class LengthModel : public ProbModel {

protected:

    std::vector< PyDist< WordId,PySparseIndex<WordId> > > sepPhrases_;
    std::vector< DirichletDist<int> > sepFor_, sepTerm_;
    std::vector< std::vector<Prob> > sepType_;
    std::vector< Prob > sepFallbacks_, sepSplits_;
    std::vector< int > phraseIdxs_, hist;
    Prob sentPen_;

public:

    LengthModel() : ProbModel() {
		
	}
	
	void setMaxLen(int maxLen) {
        sepPhrases_ = std::vector< PyDist< WordId,PySparseIndex<WordId> > >(maxLen*2, PyDist< WordId,PySparseIndex<WordId> >(1.0,0.95));
        sepFor_ = std::vector< DirichletDist<int> >(maxLen*2,DirichletDist<int>(1.0,2)) ;
        sepTerm_ = std::vector< DirichletDist<int> >(maxLen*2,DirichletDist<int>(1.0,2)); 
        sepType_ = std::vector< std::vector<Prob> >(maxLen*2,std::vector<Prob>(3,0));
        sepFallbacks_ = std::vector<Prob>(maxLen*2,0);
		sepSplits_ = std::vector<Prob>(maxLen*2);
        for(int i = 0; i < maxLen*2; i++) {
            sepPhrases_[i].setRecursive(false);
            sepSplits_[i] = -1*log(std::max(i,1));
        }
        sentPen_ = -1*log(maxLen*2);
    }
    
    Prob calcSentProb(const Span & mySpan) const { return sentPen_; }
    bool isHierarchical() { return true; }

    Prob calcGenProb(WordId jId, const Span & mySpan) const {
        int idx = mySpan.length()-1;
        return (idx || rememberNull_) ? log(sepPhrases_[idx].getProb(jId,0)) : NEG_INFINITY;
    }
    
    Prob calcBaseProb(const Span & mySpan, Prob baseMeas) const {
        int idx = mySpan.length()-1;
        return sepFallbacks_[idx]+sepType_[idx][TYPE_TERM]+baseMeas;
    }

    Prob calcTreeProb(const Span & mySpan, Prob myProb, const Span & yourSpan, Prob yourProb, int type) const {
        int idx = mySpan.length()+yourSpan.length()-1;
#ifdef DEBUG_ON
        if(idx == 0)
            throw std::runtime_error("illegal value in calcTreeProb");
#endif
        return sepFallbacks_[idx]+sepSplits_[idx]+sepType_[idx][type]+myProb+yourProb;
    }
    inline int saveIdx(WordId jId, int index) {
        while((int)phraseIdxs_.size() <= jId) 
            phraseIdxs_.resize(jId+1);
        return (phraseIdxs_[jId] = index);
    }
    void addGen(WordId jId, const Span & mySpan, Prob prob, int* tCounts, WordString & jIds) {
        int idx = saveIdx(jId,mySpan.length()-1);
        // std::cerr << "addGen("<<jId<<") --> "<<idx<<std::endl;
        sepPhrases_[saveIdx(jId,mySpan.length()-1)].addExisting(jId);
        addAverageDerivation(jId,sepPhrases_[idx].getTotal(jId),prob);
    }

    void addBase(WordId jId, const Span & mySpan, Prob prob, int* tCounts, WordString & jIds) {
        int idx = saveIdx(jId,mySpan.length()-1);
        // std::cerr << "addBase("<<jId<<") --> "<<idx<<std::endl;
        sepPhrases_[idx].addNew(jId,-1,-1,TYPE_TERM);
        addAverageDerivation(jId,sepPhrases_[idx].getTotal(jId),prob);
        addType(TYPE_TERM,idx);
    }

    void addTree(WordId jId, WordId lId, WordId rId, 
                const Span & jSpan, const Span & lSpan, const Span & rSpan, 
                int type, Prob prob, int* tCounts, WordString & jIds) {
        int idx = saveIdx(jId,jSpan.length()-1);
        // std::cerr << "addTree("<<jId<<") --> "<<idx<<std::endl;
        sepPhrases_[idx].addNew(jId, lId, rId, type);
        addAverageDerivation(jId,sepPhrases_[idx].getTotal(jId),prob);
        addType(type,idx);
    }
    
    void removeSentence(WordId head, WordString & jIds, int* counts) {
        if(head == -1)
            return;
        // std::cerr << "removeSentence: head="<<head<<","<<phraseIdxs_[head]<<std::endl;
        int idx = phraseIdxs_[head];
        // std::cerr << "removeSentence("<<head<<") --> "<<idx<<std::endl;
        PyDist< WordId,PySparseIndex<WordId> > & dist = sepPhrases_[idx];
        dist.remove(head);
        if(dist.isRemovedTable()) {
            const PyTable<WordId> & table = dist.getLastTable();
            removeType(table.type,idx);
            removeSentence(table.right,jIds,counts);
            removeSentence(table.left,jIds,counts);
        }
    }
    
    void initialize(const WordString & e, const WordString & f, 
            ParseChart & chart, WordString & jIds, int* tCounts) {
        int len = e.length()+f.length();
        for(int i = 1; i < len; i++) {
            sepFallbacks_[i] = log(sepPhrases_[i].getFallbackProb());
            sepType_[i][0] = log(sepTerm_[i].getProb(0));
            sepType_[i][1] = log(sepTerm_[i].getProb(1)*sepFor_[i].getProb(0));
            sepType_[i][2] = log(sepTerm_[i].getProb(1)*sepFor_[i].getProb(1));
        }
    }
    
    void sampleParameters(Prob defStren, Prob defDisc) {
        for(int i = 0; i < (int)sepPhrases_.size(); i++) {
            sepPhrases_[i].sampleParameters();
            if(defStren >= 0) sepPhrases_[i].setStrength(defStren);
            if(defDisc >= 0) sepPhrases_[i].setDiscount(defDisc);
        }
    }

    virtual void printStats(std::ostream &out) const {
        out << " s =";
        for(int i = 0; i < (int)sepPhrases_.size(); i++) 
            out<<" "<<sepPhrases_[i].getStrength();
        out << std::endl << " d =";
        for(int i = 0; i < (int)sepPhrases_.size(); i++)
            out<<" "<<sepPhrases_[i].getDiscount();
        out << std::endl << " t =";
        for(int i = 0; i < (int)sepPhrases_.size(); i++)
            out<<" "<<exp(sepType_[i][0])<<"/"<<exp(sepType_[i][1])<<"/"<<exp(sepType_[i][2]);
        out << std::endl;
    }


    void calcPhraseTable(const PairWordMap & jPhrases, std::vector<Prob> & eProbs, std::vector<Prob> & fProbs, std::vector<Prob> & jProbs, std::vector<Prob> & dProbs) {
        Prob myProb;
        for(PairWordMap::const_iterator it = jPhrases.begin(); it != jPhrases.end(); it++) {
            int idx = phraseIdxs_[it->second];
            if(idx || rememberNull_) {
                myProb = sepPhrases_[idx].getProb(it->second,0);
                if(myProb != 0.0) {
                    if((int)jProbs.size() <= it->second) jProbs.resize(it->second+1,0);
                    jProbs[it->second] = myProb;
                    if((int)eProbs.size() <= it->first.first) eProbs.resize(it->first.first+1,0);
                    eProbs[it->first.first] += myProb;
                    if((int)fProbs.size() <= it->first.second) fProbs.resize(it->first.second+1,0);
                    fProbs[it->first.second] += myProb;
                }
            }
        }
        dProbs = derivations_;
    }

    inline void addType(WordId id, int idx) {
        if(id == 0)
            sepTerm_[idx].add(0);
        else {
            sepTerm_[idx].add(1);
            sepFor_[idx].add(id-1);
        }
    }
    inline void removeType(WordId id, int idx) {
        if(id == 0)
            sepTerm_[idx].remove(0);
        else {
            sepTerm_[idx].remove(1);
            sepFor_[idx].remove(id-1);
        }
    }

};

}

#endif
