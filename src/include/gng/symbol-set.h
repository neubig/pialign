#ifndef GNG_SYMBOL_SET_H__
#define GNG_SYMBOL_SET_H__

#include "string.h"
#include <tr1/unordered_map>
#include <vector>

namespace gng {

template < class Key, class T, class Hash = std::tr1::hash<Key> >
class SymbolSet {

protected:
    typedef std::tr1::unordered_map< Key, T, Hash > Map;
    typedef std::vector< Key > Vocab;
    
    Map map_;
    Vocab vocab_;

public:
    SymbolSet() { }

    const Key & getSymbol(T id) const {
        return vocab_[id];
    }
    T getId(const Key & sym, bool add = false) {
        typename Map::const_iterator it = map_.find(sym);
        if(it != map_.end())
            return it->second;
        else if(add) {
            T id = vocab_.size();
            map_.insert(std::pair<Key,T>(sym,id));
            vocab_.push_back(sym);
            return id;
        }
        return -1;
    }
    size_t size() const { return vocab_.size(); }

};

}

#endif
