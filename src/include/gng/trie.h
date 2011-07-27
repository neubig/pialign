#ifndef GNG_TRIE_H__
#define GNG_TRIE_H__

// This file implements Patricia tries. They take a string of "Key"s
// and hold "Values". Values are assumed to be numeric, and negative
// values are not allowed.

#include <utility>
#include <vector>
#include <stdexcept>
#include <cstdlib>
#include <cstring>
#include <iostream>

namespace gng {

// a node of a patricia trie
template < class Key, class Value >
class TrieNode {
protected:
    typedef std::vector< std::pair< Key, TrieNode* > > ChildVec;
    Key* str_;
    int len_;
    Value val_;
    ChildVec children_;

public:
    
    // create a node with sorted stirng str of length len, and value val
    // child is added when splitting an existing node
    TrieNode(const Key* str, int len, const Value & val, const ChildVec & children = ChildVec()) : len_(len), val_(val), children_(children) {
        int size = len*sizeof(Key);
        // allocate memory
        str_ = (Key*)malloc(len*size); 
        memcpy(str_,str,len*size);
    }
    ~TrieNode() {
        free(str_);
        for(typename ChildVec::iterator it = children_.begin(); it != children_.end(); it++) {
            if(it->second) delete it->second;
        }
    }

    // find a child with binary search, add if it doesn't exists
    std::pair<Key,TrieNode*>* addChild(Key k) {
        typename ChildVec::iterator b = children_.begin(), e = children_.end(), m;
        Key comp;
        while(b != e) {
            m = b+(e-b)/2;
            comp = m->first;
            if(k > comp) b = m+1;
            else if(k < comp) e = m;
            else { return &(*m); }
        }
        b = children_.insert(b,std::pair<Key,TrieNode*>(k,0));
        return &(*b);
    }

    // find a child position with binary search, return the place
    //  it should be inserted otherwise
    typename ChildVec::iterator findChildPos(Key k) {
        typename ChildVec::iterator b = children_.begin(), e = children_.end(), m = b+(e-b)/2;
        Key comp;
        while(b != e) {
            m = b+(e-b)/2;
            comp = m->first;
            if(k > comp) b = m+1;
            else if(k < comp) e = m;
            else { return m; }
        }
        return children_.end();
    }
    
    // find a child with binary search, return null if it doesn't exist.
    // double implementation (with findChildPos) is a Bad Thing, but probably 
    // better for efficiency
    const std::pair<Key,TrieNode*>* findChild(Key k) const {
        typename ChildVec::const_iterator b = children_.begin(), e = children_.end(), m;
        Key comp;
        while(b != e) {
            m = b+(e-b)/2;
            comp = m->first;
            if(k > comp) b = m+1;
            else if(k < comp) e = m;
            else { return &(*m); }
        }
        return 0;
    }

    // insert a key
    bool insert(const Key* newStr, int newLen, Value newVal) {
        // follow the mutual strings until we can no longer
        int i, stop = std::min(len_,newLen);
        for(i = 0; i < stop && newStr[i] == str_[i]; i++);
        // if we've not reached the end of the current string, split it into
        //  a new node
        if(i < len_) {
            std::vector< std::pair< Key, TrieNode* > > oldChildren = children_;
            children_.clear();
            std::vector< std::pair< Key, TrieNode* > > newChildren(1);
            newChildren[0].first = str_[i];
            std::pair<Key,TrieNode*> * child = addChild(str_[i]);
            child->second = new TrieNode(str_+i+1,len_-i-1,val_,oldChildren);
            // intentionally do not free the string for speed
            val_ = -1; len_ = i;
        }
        // if we've reached the end of the new string, add the value
        if(i == newLen) {
            if(val_ >= 0) return false;
            else { val_ = newVal; return true; }
        }
        // otherwise, insert it into a child
        else {
            std::pair<Key,TrieNode*>* child = addChild(newStr[i]);
            if(!child->second) {
                child->second = new TrieNode(newStr+i+1,newLen-i-1,newVal);
                return true;
            } else {
                return child->second->insert(newStr+i+1,newLen-i-1,newVal);
            }
        }
    }

    // delete a key from the trie, and return whether or not the node
    //  existed
    bool erase(const Key* newStr, int newLen) {
        // follow the mutual strings until we can no longer
        int i, stop = std::min(len_,newLen);
        for(i = 0; i < stop && newStr[i] == str_[i]; i++);
        // if we've stopped halfway in between the current string, not found
        if(i != len_) return false;
        // if we've not reached the end, descend into the next node
        bool ret = false;
        if(i == newLen) {
            ret = (val_ >= 0);
            val_ = -1;
        } else {
            // search for a child and return false if not found
            typename ChildVec::iterator m = findChildPos(newStr[i]);
            if(m == children_.end()) return false;
            ret = m->second->erase(newStr+i+1,newLen-i-1);
            // if the child node is dead, remove it
            if(m->second->val_ < 0 && m->second->children_.size() == 0) {
                delete m->second;
                children_.erase(m);
            }
        }
        return ret;
    }

    Value findValue(const Key* newStr, int newLen) const {
        // follow the mutual strings until we can no longer
        int i, stop = std::min(len_,newLen);
        for(i = 0; i < stop && newStr[i] == str_[i]; i++);
        // if we've stopped halfway in between the current string, bad
        if(i < len_) return -1;
        // if we've reached the end of the string we're searching return
        if(i == newLen) return val_;
        // search for a child
        const std::pair<Key,TrieNode*> * child = findChild(newStr[i]);
        // if we found one, continue
        if(child) return child->second->findValue(newStr+i+1,newLen-i-1);
        return -1;
    }

    void commonPrefix(const Key* newStr, int newLen, int prevLen, std::vector<std::pair<int,Value> >& ret) const {
        // follow the mutual strings until we can no longer
        int i, stop = std::min(len_,newLen);
        for(i = 0; i < stop && newStr[i] == str_[i]; i++);
        // if we've stopped halfway in between the current string, bad
        if(i < len_) return;
        // if we've hit a value, add it to the array
        prevLen += len_;
        if(val_ >= 0) ret.push_back(std::pair<int,Value>(prevLen,val_));
        // if we've reached the end of the string we're searching return
        if(i == newLen) return;
        // search for a child
        const std::pair<Key,TrieNode*> * child = findChild(newStr[i]);
        // if we found one, continue
        if(child) child->second->commonPrefix(newStr+i+1,newLen-i-1,prevLen+1,ret);
    }

};

template < class Key, class Value >
class Trie {

private:

    TrieNode<Key,Value> root_;
    int size_;

public:

    Trie<Key,Value>() : root_(0, 0, -1), size_(0) { }

    // insert a value into a trie
    bool insert(const Key* newStr, int newLen, Value newVal) {
        if(newVal < 0) throw std::runtime_error("Attempt to insert invalid negative value into a trie");
        bool added = root_.insert(newStr,newLen,newVal);
        if(added) size_++;
        return added;
    }
    
    // insert a value into a trie
    bool erase(const Key* newStr, int newLen) {
        bool erased = root_.erase(newStr,newLen);
        if(erased) size_--;
        return erased;
    }
    
    // find a value from the trie
    Value findValue(const Key* newStr, int newLen) const {
        return root_.findValue(newStr,newLen);
    }

    // find all the common prefixes in the order length,value
    std::vector<std::pair<int,Value> > commonPrefix(const Key* newStr, int newLen) const {
        std::vector<std::pair<int,Value> > ret;
        root_.commonPrefix(newStr,newLen,0,ret);
        return ret;
    }

    int size() const { return size_; }
    int maxId() const { return size_; }

};

}

#endif
