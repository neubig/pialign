#include "pialign/pialign.h"
#include "pialign/model-flat.h"
#include "pialign/model-hier.h"
#include "pialign/model-length.h"

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
    WordString ws(l); 
    for(int i = 0; i < l; i++) ws[i] = i;
    return ws;
}

SpanNode * getSpanNode() {
    // make bottom words (first two are a pair)
    SpanNode * w1 = new SpanNode(Span(0,1,0,1));
    w1->type = TYPE_BASE; w1->add = false; w1->prob = -5.0;
    SpanNode * w2 = new SpanNode(Span(1,2,1,1));
    w2->type = TYPE_BASE; w2->add = false; w2->prob = -5.0;
    SpanNode * w3 = new SpanNode(Span(2,3,2,3));
    w3->type = TYPE_BASE; w3->prob = -5.0;
    SpanNode * w4 = new SpanNode(Span(3,4,1,2));
    w4->type = TYPE_BASE; w4->prob = -5.0;

    // make middle words
    SpanNode * m1 = new SpanNode(Span(0,2,0,1));
    m1->type = TYPE_BASE; m1->left = w1; m1->right = w2; m1->prob = -6.0;
    SpanNode * m2 = new SpanNode(Span(2,4,1,3));
    m2->type = TYPE_INV; m2->left = w3; m2->right = w4; m2->prob = -6.0;

    // make top words
    SpanNode * t1 = new SpanNode(Span(0,4,0,3));
    t1->type = TYPE_REG; t1->left = m1; t1->right = m2; t1->prob = -7.0;
    
    return t1;
}


// test the flat model
int testModel(ProbModel & mod, int phraseCount) {
    // make the flat model and auxiliary data structures
    mod.setMaxLen(4);
    WordString e = makeWS(4), f = makeWS(3);
    StringWordMap ep, fp;
    PairWordMap jp;

    // test adding a sample
    SpanNode * head = getSpanNode();
    mod.addSentence(e,f,head,ep,fp,jp);
    if(jp.size() != phraseCount) { 
        cerr << jp.size() << " != " << phraseCount << endl; 
        return er("testModel","unexpected number of phrases"); 
    }
    // test adding from generative dist also
    head->left->type = TYPE_GEN;
    mod.addSentence(e,f,head,ep,fp,jp);

    // test deleting the sample
    mod.removeSentence(head);
    // and delete the non-generative one
    head->left->type = TYPE_BASE;
    mod.removeSentence(head);
    if(!mod.checkEmpty())
        return er("testModel", "model was not empty");

    // clean up and return
    delete head;
    return er("testModel",0);
}

int testFlatModel() {
    FlatModel mod;
    cerr << "FlatModel: ";
    return testModel(mod, 3);
}
int testHierModel() {
    HierModel mod;
    cerr << "HierModel: ";
    return testModel(mod, 5);
}
int testLengthModel() {
    LengthModel mod;
    cerr << "LengthModel: ";
    return testModel(mod, 5);
}

int main(int argc, const char** argv) {
    int passed = 0, total = 0;

    passed += testFlatModel(); total++;
    passed += testHierModel(); total++;
    passed += testLengthModel(); total++;

    cerr << passed << "/"<<total<<" passed\n";

}
