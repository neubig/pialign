#ifndef PIALIGN_PORT_PORT_H_
#define PIALIGN_PORT_PORT_H_

// As of OS X 10.9, it looks like C++ TR1 headers are removed from the
// search paths. Instead, we can include C++11 headers.
#if defined(__APPLE__)
#include <AvailabilityMacros.h>
#endif

#include <functional>
#include <unordered_map>
#include <unordered_set>

#endif // PIALIGN_PORT_PORT_H_
