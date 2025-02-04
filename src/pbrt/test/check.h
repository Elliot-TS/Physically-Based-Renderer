#pragma once
#include <cassert>
#define CHECK(x) assert(x)
#define DCHECK(x) (CHECK(x))
#define DCHECK_NE(a, b) (CHECK(a != b))

#define LOG(x) std::cout << x << std::endl

// Log Move Semantics
#ifndef LOG_MOVE
#define LOG_MOVE(x)
#endif
#define CREATED(x) LOG_MOVE(x << " Created")
#define MOVED(x) LOG_MOVE(x << " Moved")
#define DELETED(x) LOG_MOVE(x << " Deleted")
#define COPIED(x) LOG_MOVE(x << " Copied")

