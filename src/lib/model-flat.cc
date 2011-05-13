
#include "pialign/model-flat.h"

using namespace pialign;

void FlatModel::addSentence(const WordString & e, const WordString & f, SpanNode* node, StringWordMap & ePhrases, StringWordMap & fPhrases, PairWordMap & pairs) {
    if(!node || !node->add) return;
    // if we're on a non-terminal, add the type and recurse
    if(node->type == TYPE_REG || node->type == TYPE_INV) {
        addType(node->type);
        addSentence(e,f,node->left, ePhrases,fPhrases,pairs);
        addSentence(e,f,node->right,ePhrases,fPhrases,pairs);
    } else {
        // add the terminal symbol
        addType(TYPE_TERM);
        // find the phrase ID
        const Span & mySpan = node->span;
        int s=mySpan.es,t=mySpan.ee,u=mySpan.fs,v=mySpan.fe;
        WordId eId = ePhrases.getId(e.substr(s,t-s),true),
            fId = fPhrases.getId(f.substr(u,v-u),true);
        // set the phrase and the pair
        node->phraseid = pairs.getId(std::pair<WordId,WordId>(eId,fId),true);
        if(node->type == TYPE_BASE)
            phrases_.addNew(node->phraseid,-1,-1,-1);
        else
            phrases_.addExisting(node->phraseid);
        addAverageDerivation(node->phraseid,phrases_.getTotal(node->phraseid),node->prob);
    }
}

void FlatModel::removeSentence(const SpanNode* node) {
    if(!node || !node->add) return;
    if(node->type == TYPE_REG || node->type == TYPE_INV) {
        removeType(node->type);
        // std::cerr << " T("<<node->type<<")";
        removeSentence(node->left);
        removeSentence(node->right);
    }
    else {
        removeType(TYPE_TERM);
        phrases_.remove(node->phraseid);
        // std::cerr << " P(" << node->span << ") T("<<TYPE_TERM<<")";
    }
}
