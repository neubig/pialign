
#include "model-length.h"

using namespace std;
using namespace std::tr1;
using namespace pialign;
using namespace gng;


void LengthModel::setMaxLen(int maxLen) {
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


void LengthModel::addSentence(const WordString & e, const WordString & f, SpanNode* node, StringWordMap & ePhrases, StringWordMap & fPhrases, PairWordMap & pairs) {
    if(!node || !node->add) return;
    // get the phrase IDs
    const Span & mySpan = node->span;
    int s=mySpan.es,t=mySpan.ee,u=mySpan.fs,v=mySpan.fe;
    WordId eId = ePhrases.getId(e.substr(s,t-s),true),
        fId = fPhrases.getId(f.substr(u,v-u),true);
    node->phraseid = pairs.getId(std::pair<WordId,WordId>(eId,fId),true);
    int idx = saveIdx(node->phraseid,node->span.length()-1);
    int toAdd = node->type;
    // handle either non-terminals or terminals
    if(toAdd == TYPE_REG || toAdd == TYPE_INV) {
        addSentence(e,f,node->left, ePhrases,fPhrases,pairs);
        addSentence(e,f,node->right,ePhrases,fPhrases,pairs);
    } else
        toAdd = TYPE_TERM;
    addType(toAdd,idx);
    // find the left and right nodes
    WordId lId = (node->left?node->left->phraseid:-1),
            rId = (node->right?node->right->phraseid:-1);
    // add the appropriate values for the derivation
    if(node->type == TYPE_GEN)
        sepPhrases_[idx].addExisting(node->phraseid);
    else
        sepPhrases_[idx].addNew(node->phraseid,lId,rId,toAdd);
    addAverageDerivation(node->phraseid,sepPhrases_[idx].getTotal(node->phraseid),node->prob);
}


void LengthModel::removePhrasePair(WordId jId) {
    if(jId < 0) return;
#ifdef DEBUG_ON
    if(jId >= phraseIdxs_.size()) {
        std::cerr << jId << " >= " << phraseIdxs_.size() << std::endl;
        throw std::runtime_error("Overflown phraseIdx in model-length.h");
    }
#endif
    int idx = phraseIdxs_[jId];
    PyDist< WordId,PySparseIndex<WordId> > & dist = sepPhrases_[idx];
    dist.remove(jId);
    if(dist.isRemovedTable()) {
        const PyTable<WordId> & table = dist.getLastTable();
        removeType(table.type,idx);
        removePhrasePair(table.right);
        removePhrasePair(table.left);
    }
}


void LengthModel::initialize(const WordString & e, const WordString & f, 
        ParseChart & chart, WordString & jIds, int* tCounts) {
    int len = e.length()+f.length();
    for(int i = 1; i < len; i++) {
        sepFallbacks_[i] = log(sepPhrases_[i].getFallbackProb());
        sepType_[i][0] = log(sepTerm_[i].getProb(0));
        sepType_[i][1] = log(sepTerm_[i].getProb(1)*sepFor_[i].getProb(0));
        sepType_[i][2] = log(sepTerm_[i].getProb(1)*sepFor_[i].getProb(1));
    }
}


void LengthModel::printStats(std::ostream &out) const {
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


void LengthModel::calcPhraseTable(const PairWordMap & jPhrases, std::vector<Prob> & eProbs, std::vector<Prob> & fProbs, std::vector<Prob> & jProbs, std::vector<Prob> & dProbs) {
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
