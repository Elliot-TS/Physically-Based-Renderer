#pragma once
#include <cassert>
#define CHECK(x)        assert(x)
#define CHECK_NE(a, b)  (CHECK(a != b))
#define CHECK_EQ(a, b)  (CHECK(a == b))
#define DCHECK(x)       (CHECK(x))
#define DCHECK_NE(a, b) (DCHECK(a != b))
#define DCHECK_GT(a, b) (DCHECK(a > b))
#define DCHECK_LT(a, b) (DCHECK(a < b))
#define DCHECK_GE(a, b) (DCHECK(a >= b))
#define DCHECK_LE(a, b) (DCHECK(a <= b))
#define LOG_FATAL(s)    std::cerr << s << std::endl

#define LOG(x)   std::cout << x << std::endl
#define ERROR(x) std::cerr << x << std::endl

#define ASSERT(condition, message) \
  [&]() {                          \
    bool c = condition;            \
    if (!c) ERROR(message);        \
    return c;                      \
  }()
// #define ASSERT(condition, message) condition

// Log Move Semantics
#ifndef LOG_MOVE
  #define LOG_MOVE(x)
#endif
#define CREATED(x) LOG_MOVE(x << " Created")
#define MOVED(x)   LOG_MOVE(x << " Moved")
#define DELETED(x) LOG_MOVE(x << " Deleted")
#define COPIED(x)  LOG_MOVE(x << " Copied")
