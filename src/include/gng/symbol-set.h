#ifndef GNG_SYMBOL_SET_H__
#define GNG_SYMBOL_SET_H__

#include "string.h"
#include <tr1/unordered_map>
#include <vector>

namespace gng {

template < class Key, class T, class Hash = std::tr1::hash<Key> >
class SymbolSet {

public:

    typedef std::tr1::unordered_map< Key, T, Hash > Map;
    typedef std::vector< Key > Vocab;
    typedef typename Map::iterator iterator;
    typedef typename Map::const_iterator const_iterator;

protected:
    
    Map map_;
    Vocab vocab_;
    std::vector<T> nextIds_;

public:
    SymbolSet() { }

    unsigned numElements() { return vocab_.size() - nextIds_.size(); }

    const Key & getSymbol(T id) const {
        return vocab_[id];
    }
    T getId(const Key & sym, bool add = false) {
        typename Map::const_iterator it = map_.find(sym);
        if(it != map_.end())
            return it->second;
        else if(add) {
            T id;
            if(nextIds_.size()) { 
                id = nextIds_[nextIds_.size()-1];
                nextIds_.pop_back();
                vocab_[id] = sym;
            }
            else {
                id = vocab_.size();
                vocab_.push_back(sym);
            }
            map_.insert(std::pair<Key,T>(sym,id));
            return id;
        }
        return -1;
    }
    T getId(const Key & sym) const {
        typename Map::const_iterator it = map_.find(sym);
        if(it != map_.end())
            return it->second;
        return -1;
    }
    size_t size() const { return vocab_.size(); }


    void removeElements(const std::vector<T> & vec) {
        for(unsigned i = 0; i < vec.size(); i++) {
            map_.erase(vocab_[vec[i]]);
            vocab_[vec[i]] = Key();
            nextIds_.push_back(vec[i]);
        }
    }
    
    typename Map::iterator begin() { return map_.begin(); }
    typename Map::iterator end() { return map_.end(); }
    typename Map::const_iterator begin() const { return map_.begin(); }
    typename Map::const_iterator end() const { return map_.end(); }

};

}

#endif
