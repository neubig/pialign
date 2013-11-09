#include "gng/trie.h"
#include <iostream>

using namespace gng;
using namespace std;

int er(const char* func, const char* er) {
    if(er) {
        cerr << func << ": FAILED with "<<er << endl;
        return 0;
    }
    cerr << func << ": OK" << endl;
    return 1;
}

int testInsert() {

    Trie<char,int> trie;

    // insert new
    if(!trie.insert("abracadabra",11,1)) { return er("testInsert", "abracadabra"); }
    if(!trie.insert("abs",3,2)) { return er("testInsert", "abs"); }
    if(!trie.insert("babies",6,3)) { return er("testInsert", "babies"); }
    if(!trie.insert("baby",4,4)) { return er("testInsert", "baby"); }
    if(!trie.insert("babino",6,5)) { return er("testInsert", "babino"); }
    if(!trie.insert("babel",5,6)) { return er("testInsert", "babel"); }

    // can't insert again
    if(trie.insert("abracadabra",11,1)) { return er("testInsert", "abracadabra reinsert"); }
    if(trie.insert("abs",3,2)) { return er("testInsert", "abs reinsert"); }
    if(trie.insert("babies",6,3)) { return er("testInsert", "babies reinsert"); }
    if(trie.insert("baby",4,4)) { return er("testInsert", "baby reinsert"); }
    if(trie.insert("babino",6,5)) { return er("testInsert", "babino reinsert"); }
    if(trie.insert("babel",5,6)) { return er("testInsert", "babel reinsert"); }

    return er("testInsert", 0);

}

int testFind() {

    Trie<char,int> trie;

    // insert new
    if(!trie.insert("abracadabra",11,1)) { return er("testFind", "insert abracadabra"); }
    if(!trie.insert("abs",3,2)) { return er("testFind", "insert abs"); }
    if(!trie.insert("babies",6,3)) { return er("testFind", "insert babies"); }
    if(!trie.insert("baby",4,4)) { return er("testFind", "insert baby"); }
    if(!trie.insert("babino",6,5)) { return er("testFind", "insert babino"); }
    if(!trie.insert("babel",5,6)) { return er("testFind", "insert babel"); }

    if(trie.findValue("abracadabra",11) != 1) { return er("testFind", "find abracadabra"); }
    if(trie.findValue("abs",3) != 2) { return er("testFind", "find abs"); }
    if(trie.findValue("babies",6) != 3) { return er("testFind", "find babies"); }
    if(trie.findValue("baby",4) != 4) { return er("testFind", "find baby"); }
    if(trie.findValue("babino",6) != 5) { return er("testFind", "find babino"); }
    if(trie.findValue("babel",5) != 6) { return er("testFind", "find babel"); }
    if(trie.findValue("bab",3) != -1) { return er("testFind", "find bab"); }
    if(trie.findValue("babys",5) != -1) { return er("testFind", "find babys"); }

    return er("testFind", 0);

}

int testErase() {

    Trie<char,int> trie;

    // insert new
    if(!trie.insert("abracadabra",11,1)) { return er("testErase", "insert abracadabra"); }
    if(!trie.insert("abs",3,2)) { return er("testErase", "insert abs"); }
    if(!trie.insert("babies",6,3)) { return er("testErase", "insert babies"); }
    if(!trie.insert("baby",4,4)) { return er("testErase", "insert baby"); }
    if(!trie.insert("babino",6,5)) { return er("testErase", "insert babino"); }
    if(!trie.insert("babel",5,6)) { return er("testErase", "insert babel"); }

    if(trie.size() != 6) { return er("testErase", "size after insert"); }

    if(!trie.erase("abs",3)) { return er("testErase", "erase abs"); }
    if(!trie.erase("babel",5)) { return er("testErase", "erase babel"); }
    if(trie.erase("babel",5)) { return er("testErase", "reerase babel"); }
    if(trie.erase("chorus",6)) { return er("testErase", "erase chorus"); }

    if(trie.size() != 4) { return er("testErase", "size after erase"); }

    if(trie.findValue("abracadabra",11) != 1) { return er("testErase", "find abracadabra"); }
    if(trie.findValue("abs",3) != -1) { return er("testErase", "find abs"); }
    if(trie.findValue("babies",6) != 3) { return er("testErase", "find babies"); }
    if(trie.findValue("baby",4) != 4) { return er("testErase", "find baby"); }
    if(trie.findValue("babino",6) != 5) { return er("testErase", "find babino"); }
    if(trie.findValue("babel",5) != -1) { return er("testErase", "find babel"); }
    if(trie.findValue("bab",3) != -1) { return er("testErase", "find bab"); }
    if(trie.findValue("babys",5) != -1) { return er("testErase", "find babys"); }

    return er("testErase", 0);

}


int main(int argc, const char** argv) {
    int passed = 0, total = 0;

    passed += testInsert(); total++;
    passed += testFind(); total++;
    passed += testErase(); total++;

    cerr << passed << "/"<<total<<" passed\n";

}
