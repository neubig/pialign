
#include "gng/counter.h"
#include "pialign/base-phrasecooc.h"
#include "pialign/esa.hxx"
#include <map>

using namespace std;
using namespace pialign;
using namespace gng;

void BasePhraseCooc::substringMatrix(Corpus & corp, const WordSymbolSet & vocab, int boost, vector<WordString> & subs, vector< set<int> > & sents, Prob coocDisc) {
    // make the suffix tree
    vector<WordId> & T = corp.getData();
    vector<int> SA(T.size());
    vector<int> L (T.size());
    vector<int> R (T.size());
    vector<int> D (T.size());
    int nodeNum = 0;
    vector<int> sids = corp.getSentIds();
    int numSents = sids[sids.size()-1];
    if (esaxx(T.begin(), SA.begin(), 
          L.begin(), R.begin(), D.begin(), 
          (int)T.size(), (int)vocab.size()+boost+numSents+1, nodeNum) == -1){
        throw std::runtime_error("could not make suffix tree");
    }
    // add every string that appears at least twice, and its sentences
    for (int i = 0; i < nodeNum; ++i) {
        int beg = SA[L[i]];
        if(D[i] > 0) {
            set<int> mySents;
            // check to make sure the left size is not consistent
            bool leftDiff = false;
            int type = (SA[L[i]]==0||sids[SA[L[i]]]!=sids[SA[L[i]]-1])?-1:T[SA[L[i]]-1];
            for(int j = L[i]; j < R[i]; j++) {
                int myType = (SA[j]==0||sids[SA[j]]!=sids[SA[j]-1])?-1:T[SA[j]-1];
                leftDiff = leftDiff || myType != type || myType == -1;
                if(sids[SA[j]] == sids[SA[j]+D[i]-1] && corp[sids[SA[j]]].length() != 0)
                    mySents.insert(sids[SA[j]]);
            }
            if(leftDiff && mySents.size() > coocDisc) {
                subs.push_back(WordString(&T[beg],D[i]));
                sents.push_back(mySents);
            }
        }
    }
    cerr << "T.size="<<T.size()<<", vocab.size()="<<vocab.size()<<", nodes="<<nodeNum<<", words="<<subs.size()<<endl;
}

void BasePhraseCooc::trainCooc(Corpus & es, const WordSymbolSet & eVocab, Corpus & fs, const WordSymbolSet & fVocab, Prob coocDisc, Prob coocCut) {
    vector< WordString > eStrs, fStrs;
    vector< set<int> > eSents, fSents;
    // calculate the substring sparse matrix in e
    cerr << "Calculating substrings for e"<<endl;
    substringMatrix(es,eVocab,0,eStrs,eSents, coocDisc);
    // calculate the substring sparse matrix in f
    cerr << "Calculating substrings for f"<<endl;
    substringMatrix(fs,fVocab,eVocab.size(),fStrs,fSents, coocDisc);
    // create the inverted matrix for f
    cerr << "Creating inverted matrix for f"<<endl;
    vector< vector<int> > fWords(fs.size());
    for(unsigned i = 0; i < fSents.size(); i++)
        for(set<int>::const_iterator it = fSents[i].begin(); it != fSents[i].end(); it++) {
            fWords[*it].push_back(i);
        }
    // take the matrix product and calculate the Cooc coefficient for all variables
    cerr << "Calculating co-occurence"<<endl;
    Prob sum = 0;
    vector<WordId> eIds(eSents.size(), -1), fIds(fSents.size(), -1);
    for(unsigned i = 0; i < eSents.size(); i++) {
        Counter<int,int> myCounts;
        for(set<int>::const_iterator eit = eSents[i].begin(); eit != eSents[i].end(); eit++) {
            for(vector<int>::const_iterator fit = fWords[*eit].begin(); fit != fWords[*eit].end(); fit++) {
                myCounts.inc(*fit);
            }
        }
        for(Counter<int,int>::const_iterator cit = myCounts.begin(); cit != myCounts.end(); cit++) {
            if(cit->second > coocDisc) {
                if(eIds[i] == -1) eIds[i] = eSymbols_.getId(eStrs[i],true);
                if(fIds[cit->first] == -1) fIds[cit->first] = fSymbols_.getId(fStrs[cit->first],true);
                int eCount = eSents[i].size(), fCount = fSents[cit->first].size(), efCount = cit->second;
                double eProb = (eCount-coocDisc)/es.size(), fProb = (fCount-coocDisc)/es.size(), efProb = (efCount-coocDisc)/es.size();
                double geomMean = efProb/eProb*efProb/fProb;
                sum += geomMean;
                jProbs_.insert(PairProbMap::value_type(pair<WordId,WordId>(eIds[i],fIds[cit->first]),geomMean));

                // // ---- print ----
                // for(unsigned j = 0; j < eStrs[i].length(); j++) {
                //     if(j != 0) cout << " ";
                //     cout << eVocab.getSymbol(eStrs[i][j]);
                // }
                // cout << "\t";
                // for(unsigned j = 0; j < fStrs[cit->first].length(); j++) {
                //     if(j != 0) cout << " ";
                //     cout <<fVocab.getSymbol(fStrs[cit->first][j]-eVocab.size());
                // }
                // cout << "\t"<<efCount<<"\t"<<eCount<<"\t"<<fCount<<endl;
            }
        }
    }
    // normalize probabilities and remove low probability values
    cerr << "Normalizing and cleaning probabilities " << endl;
    vector< pair<WordId,WordId> > remove;
    uniPen_ = 0;
    for(PairProbMap::iterator it = jProbs_.begin(); it != jProbs_.end(); it++) {
        it->second /= sum;
        if(it->second < coocCut) {
            remove.push_back(it->first);
            uniPen_ += it->second;
        }
        it->second = log(it->second);
    }
    for(unsigned i = 0; i < remove.size(); i++)
        jProbs_.erase(remove[i]);
    if(uniPen_ == 0)
        cerr << "WARNING: Base measure unigram probability is 0, try setting -cooccut higher" << endl;
    else {
        cerr << "Unigram probability "<<uniPen_<<endl;
        uniPen_ = log(uniPen_);
    }
}

void BasePhraseCooc::addBases(const WordString & e, const WordString & f, const ProbModel & mod, ParseChart & chart, SpanProbMap & baseChart) const {

    BaseUnigram::addBases(e,f,mod,chart,baseChart);
 
    // add the span probabilities
    vector< LabeledEdge > eEdges = eSymbols_.findEdges(e,e.length());
    vector< LabeledEdge > fEdges = fSymbols_.findEdges(f,f.length());
    for(int i = 0; i < (int)eEdges.size(); i++) {
        const LabeledEdge & ee = eEdges[i];
        for(int j = 0; j < (int)fEdges.size(); j++) {
            const LabeledEdge & fe = fEdges[j];
            PairProbMap::const_iterator it = 
                jProbs_.find(pair<WordId,WordId>(ee.l,fe.l));
            if(it != jProbs_.end()) { 
                // cerr << " found in jProbs (["<<ee.s<<","<<ee.e<<","<<ee.l<<"]/["<<fe.e<<","<<fe.s<<","<<fe.l<<"]) = "<<it->second<<endl;
                Span s(ee.s, ee.e, fe.s, fe.e);
                baseChart.insertProb(s, chart.addToChart(s,it->second));
            }
        }
    }
}
