#ifndef GNG_SYMBOL_TRIE_H__
#define GNG_SYMBOL_TRIE_H__

#include "string.h"
#include "gng/trie.h"
#include <vector>
#include <iostream>

namespace gng {

template < class K, class T >
class SymbolTrie {

protected:

    std::vector<T> nextKeys_;
    std::vector<std::pair<K*,int> > strs_;
    gng::Trie<K,T> map_;

public:
    SymbolTrie() { }
    ~SymbolTrie() { 
        for(unsigned i = 0; i < strs_.size(); i++)
            if(strs_[i].first)
                free(strs_[i].first);
    }

    unsigned maxId() const { return map_.size() + nextKeys_.size(); }
    unsigned numElements() const { return map_.size(); }

    T getId(const K* ptr, int len, bool add) {
        // for(int i = 0; i < len; i++) 
        //     std::cerr << ptr[i] << " ";
        // std::cerr << ", l="<<len<<", getId("<<(int)ptr<<","<<len<<","<<add<<")"<<std::endl;
        T ret = map_.findValue(ptr,len);
        if(ret < 0 && add) {
            if(nextKeys_.size()) { 
                ret = nextKeys_[nextKeys_.size()-1];
                nextKeys_.pop_back();
            }
            else {
                ret = strs_.size();
                strs_.push_back(std::pair<K*,int>(0,0));
            }
            map_.insert(ptr,len,ret);
            strs_[ret].first = (K*)malloc(len*sizeof(K)); 
            strs_[ret].second = len;
            memcpy(strs_[ret].first,ptr,len*sizeof(K)); 
        }
        return ret;
    }
    T getId(const GenericString<K> & sym, bool add = false) {
        return getId(&sym[0], sym.length(), add);
    }
    const T getId(const K* ptr, int len) const {
        return map_.findValue(ptr,len);
    }

    const K* getSymbol(T id) const { return strs_[id].first; }
    int getSymbolLength(T id) const { return strs_[id].second; }

    inline size_t size() const { return map_.size(); }
    
    void removeElements(const std::vector<T> & vec) {
        // don't do anything for now...
        for(typename std::vector<T>::const_iterator it = vec.begin(); it != vec.end(); it++) {
            nextKeys_.push_back(*it);
            map_.erase(strs_[*it].first,strs_[*it].second);
            free(strs_[*it].first); strs_[*it].first = 0;
            strs_[*it].second = 0;
        }
    }



};

}

#endif
