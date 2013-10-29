#ifndef PIALIGN_PORT_PORT_H_
#define PIALIGN_PORT_PORT_H_

// As of OS X 10.9, it looks like C++ TR1 headers are removed from the
// search paths. Instead, we can include C++11 headers.
#if defined(__APPLE__) && MAC_OS_X_VERSION_MIN_REQUIRED >= MAC_OS_X_VERSION_10_9
#include <functional>
#include <unordered_map>
#include <unordered_set>
#else // Assuming older OS X, Linux or similar platforms
#include <tr1/functional>
#include <tr1/unordered_map>
#include <tr1/unordered_set>

namespace std {
using tr1::hash;
using tr1::unordered_map;
using tr1::unordered_set;
} // namespace std
#endif

#endif // PIALIGN_PORT_PORT_H_
