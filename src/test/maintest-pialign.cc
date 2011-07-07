#include "pialign/pialign.h"
#include "pialign/model-flat.h"
#include "pialign/model-hier.h"
#include "pialign/model-length.h"
#include "pialign/base-unigram.h"

using namespace pialign;
using namespace std;

int er(const char* func, const char* err) {
    if(err) {
        cerr << func << ": FAILED with "<<err << endl;
        return 0;
    }
    cerr << func << ": OK" << endl;
    return 1;
}

WordString makeWS(int l) {
    WordId *wid = new WordId[l];
    WordString ws(wid,l); 
    for(int i = 0; i < l; i++) wid[i] = i;
    return ws;
}

SpanNode * getSpanNode() {
    // make bottom words (first two are a pair)
    SpanNode * w1 = new SpanNode(Span(0,1,0,1));
    w1->type = TYPE_BASE; w1->add = false; w1->prob = -5.0; w1->baseProb = -10.0;
    SpanNode * w2 = new SpanNode(Span(1,2,1,1));
    w2->type = TYPE_BASE; w2->add = false; w2->prob = -5.0; w2->baseProb = -10.0;
    SpanNode * w3 = new SpanNode(Span(2,3,2,3));
    w3->type = TYPE_BASE; w3->prob = -5.0; w3->baseProb = -10.0;
    SpanNode * w4 = new SpanNode(Span(3,4,1,2));
    w4->type = TYPE_BASE; w4->prob = -5.0; w4->baseProb = -10.0;

    // make middle words
    SpanNode * m1 = new SpanNode(Span(0,2,0,1));
    m1->type = TYPE_BASE; m1->left = w1; m1->right = w2; m1->prob = -6.0; m1->baseProb = -12.0;
    SpanNode * m2 = new SpanNode(Span(2,4,1,3));
    m2->type = TYPE_INV; m2->left = w3; m2->right = w4; m2->prob = -6.0;

    // make top words
    SpanNode * t1 = new SpanNode(Span(0,4,0,3));
    t1->type = TYPE_REG; t1->left = m1; t1->right = m2; t1->prob = -7.0;
    
    return t1;
}


// test various models
int testModel(ProbModel & mod, int phraseCount) {
    mod.setMaxLen(4);
    WordString e = makeWS(4), f = makeWS(3);
    BaseUnigram base;
    StringWordSet ep, fp;
    PairWordSet jp;

    // test adding a sample
    SpanNode * head = getSpanNode();
    mod.addSentence(e,f,head,ep,fp,jp,&base);
    if((int)jp.size() != phraseCount) { 
        cerr << jp.size() << " != " << phraseCount << endl; 
        return er("testModel","unexpected number of phrases"); 
    }
    // test adding from generative dist also
    head->left->type = TYPE_GEN;
    mod.addSentence(e,f,head,ep,fp,jp,&base);

    // test deleting the sample
    mod.removeSentence(head,&base);
    // and delete the non-generative one
    head->left->type = TYPE_BASE;
    mod.removeSentence(head,&base);
    if(!mod.checkEmpty())
        return er("testModel", "model was not empty");

    // clean up and return
    delete head;
    delete [] e.getPointer();
    delete [] f.getPointer();
    return er("testModel",0);
}


// test various models
int testMatch(ProbModel & mod) {
    mod.setMaxLen(4);
    WordString e = makeWS(4), f = makeWS(3);
    BaseUnigram base;
    StringWordSet ep, fp;
    PairWordSet jp;

    int size = 4;

    vector<Prob> addProbs(size,0), remProbs(size,0);
    SpanNode * head = getSpanNode();

    Prob addProb = 0, remProb = 0;

    // add the probabilities
    for(int i = 0; i < size; i++) {
        head->left->type = (i < size/2 ? TYPE_BASE : TYPE_GEN);
        // cerr << "---- adding sentence "<<i<<" ----"<<endl;
        addProbs[i] = mod.addSentence(e,f,head,ep,fp,jp,&base);
        addProb += addProbs[i];
    }
    
    // subtract the probabilities
    for(int i = size-1; i >= 0; i--) {
        head->left->type = (i < size/2 ? TYPE_BASE : TYPE_GEN);
        // cerr << "---- removing sentence "<<i<<" ----"<<endl;
        SpanNode * node = mod.removeSentence(head,&base);
        remProbs[i] = node->prob;
        delete node;
        remProb += remProbs[i];
    }

    if(!mod.checkEmpty())
        return er("testMatch", "model was not empty");

    for(int i = 0; i < size; i++) {
        if(abs(addProbs[i] - remProbs[i]) > 0.01) {
            cerr << "addProbs["<<i<<"] "<<addProbs[i]<<" != remProbs["<<i<<"] "<<remProbs[i]<<endl;
            // return er("testMatch","probabilities don't match");
        } 
    }

    if(addProb != remProb) {
        cerr << "addProb "<<addProb<<" != remProb "<<remProb<<endl;
        return er("testMatch", "probabilities did not match");
    }
    

    // clean up and return
    delete head;
    delete [] e.getPointer();
    delete [] f.getPointer();
    return er("testMatch",0);
}

int testFlatModel() {
    FlatModel mod;
    cerr << "FlatModel simple: ";
    return testModel(mod, 3);
}
int testHierModel() {
    HierModel mod;
    cerr << "HierModel simple: ";
    return testModel(mod, 5);
}
int testLengthModel() {
    LengthModel mod;
    cerr << "LengthModel simple: ";
    return testModel(mod, 5);
}

int testFlatMatch() {
    FlatModel mod;
    cerr << "FlatModel match: ";
    return testMatch(mod);
}
int testHierMatch() {
    HierModel mod;
    cerr << "HierModel match: ";
    return testMatch(mod);
}
int testLengthMatch() {
    LengthModel mod;
    cerr << "LengthModel match: ";
    return testMatch(mod);
}


int main(int argc, const char** argv) {
    int passed = 0, total = 0;

    passed += testFlatModel(); total++;
    passed += testHierModel(); total++;
    passed += testLengthModel(); total++;
    passed += testFlatMatch(); total++;
    passed += testHierMatch(); total++;
    passed += testLengthMatch(); total++;

    cerr << passed << "/"<<total<<" passed\n";

}
