
#include "pialign/model-flat.h"

using namespace pialign;
using namespace std;

Prob FlatModel::addSentence(const WordString & e, const WordString & f, SpanNode* node, StringWordSet & ePhrases, StringWordSet & fPhrases, PairWordSet & pairs, BaseMeasure* base) {
    if(!node || !node->add) return 0;
    // if we're on a non-terminal, add the type and recurse
    Prob totProb = 0;
    if(node->type == TYPE_REG || node->type == TYPE_INV) {
        totProb += addType(node->type);
        totProb += addSentence(e,f,node->right,ePhrases,fPhrases,pairs,base);
        totProb += addSentence(e,f,node->left, ePhrases,fPhrases,pairs,base);
    } else {
        // add the terminal symbol
        totProb += addType(TYPE_TERM);
        // find the phrase ID
        const Span & mySpan = node->span;
        int s=mySpan.es,t=mySpan.ee,u=mySpan.fs,v=mySpan.fe;
        WordId eId = ePhrases.getId(e.substr(s,t-s),true),
            fId = fPhrases.getId(f.substr(u,v-u),true);
        // set the phrase and the pair
        node->phraseid = pairs.getId(std::pair<WordId,WordId>(eId,fId),true);
        if(node->baseProb != 0) {
            base->add(node->span,node->phraseid,node->baseProb);
        }
        if(node->type == TYPE_BASE) {
            // std::cerr << "baseProb = "<<node->baseProb<<" + "<<log(phrases_.getFallbackProb())<<std::endl;
            totProb += node->baseProb + log(phrases_.getFallbackProb());
            phrases_.addNew(node->phraseid,-1,-1,-1);
        } else {
            // cerr << "totProb += log("<<phrases_.getProb(node->phraseid,0)<<"), "<<log(phrases_.getProb(node->phraseid,0))<<endl;
            totProb += log(phrases_.getProb(node->phraseid,0));
            phrases_.addExisting(node->phraseid);
        }
        addAverageDerivation(node->phraseid,phrases_.getTotal(node->phraseid),node->prob);
    }
    // cerr << " add: s="<<node->span<<", i="<<node->phraseid<<", t="<<node->type<<", p="<<node->prob<<", b="<<node->baseProb<<", a="<<node->add<<endl;
    return totProb;
}

SpanNode* FlatModel::removeSentence(const SpanNode* node, BaseMeasure* base) {
    if(!node || !node->add) return 0;
    SpanNode * ret = new SpanNode(Span(0,0,0,0));
    ret->phraseid = node->phraseid;
    if(node->type == TYPE_REG || node->type == TYPE_INV) {
        ret->type = node->type;
        ret->prob += removeType(node->type);
        ret->left = removeSentence(node->left,base);
        ret->right = removeSentence(node->right,base);
        ret->prob += ret->left->prob; ret->prob += ret->right->prob; 
    }
    else {
        ret->prob += removeType(TYPE_TERM);
        ret->prob += phrases_.remove(node->phraseid);
        if(phrases_.isRemovedTable()) {
            ret->prob += base->getBase(node->phraseid);
            ret->type = TYPE_BASE;
        } 
        else
            ret->type = TYPE_GEN;
    }
    // cerr << " remove: s="<<node->span<<", i="<<node->phraseid<<", t="<<node->type<<", p="<<node->prob<<", b="<<node->baseProb<<", a="<<node->add<<endl;
    return ret;
}
