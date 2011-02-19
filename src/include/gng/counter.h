#ifndef COUNTER_H__
#define COUNTER_H__

#include <tr1/unordered_map>
#include <stdexcept>
#include <string>

namespace gng {

template <class Key,class T>
class Counter : public std::tr1::unordered_map<Key,T> {


public:
    typedef typename Counter<Key,T>::iterator iterator;

    T inc(const Key & s, T weight=1) {
        iterator it = find(s);
        if(it == this->end()) {
            if(weight < 0)
                throw std::runtime_error("Negative counter");
            this->insert(std::pair<Key,T>(s,weight));
            return weight;
        } else {
            it->second += weight;
            if(it->second < 0)
                throw std::runtime_error("Negative counter");
            return it->second;
        }
    }

    T get(const Key & s) {
        iterator it = this->find(s);
        return (it == this->end()?(T)0:it->second);
    }
    
};

}

#endif
