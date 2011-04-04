#include "model-hier.h"

using namespace pialign;

void HierModel::addSentence(const WordString & e, const WordString & f, SpanNode* node, StringWordMap & ePhrases, StringWordMap & fPhrases, PairWordMap & pairs) {
    if(!node || !node->add) return;
    // get the phrase IDs
    const Span & mySpan = node->span;
    int s=mySpan.es,t=mySpan.ee,u=mySpan.fs,v=mySpan.fe;
    WordId eId = ePhrases.getId(e.substr(s,t-s),true),
        fId = fPhrases.getId(f.substr(u,v-u),true);
    node->phraseid = pairs.getId(std::pair<WordId,WordId>(eId,fId),true);
    int toAdd = node->type;
    // handle either non-terminals or terminals
    if(toAdd == TYPE_REG || toAdd == TYPE_INV) {
        addSentence(e,f,node->left,ePhrases,fPhrases,pairs);
        addSentence(e,f,node->right,ePhrases,fPhrases,pairs);
    } else
        toAdd = TYPE_TERM;
    // find the left and right nodes
    WordId lId = (node->left?node->left->phraseid:-1),
            rId = (node->right?node->right->phraseid:-1);
    // add the appropriate values for the derivation
    if(node->type == TYPE_GEN)
        phrases_.addExisting(node->phraseid);
    else {
        addType(toAdd);
        phrases_.addNew(node->phraseid,lId,rId,toAdd);
    }
    addAverageDerivation(node->phraseid,phrases_.getTotal(node->phraseid),node->prob);
}


void HierModel::removeSentence(const SpanNode* node) {
    if(!node) return;
    std::vector<int> ids = phrases_.remove(node->phraseid);
    // std::cerr << " phrases_.isRecursive? "<<phrases_.isRecursive()<<", removing "<<head<<",";
    for(int i = 0; i < (int)ids.size(); i++) {
        // std::cerr << " " << ids[i];
        removeType(ids[i]);
    }
    // std::cerr<<std::endl;
}
