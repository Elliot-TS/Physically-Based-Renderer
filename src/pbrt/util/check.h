#include <cassert>
#define CHECK(x) assert(x)
#define DCHECK(x) (CHECK(x))
#define DCHECK_NE(a, b) (CHECK(a != b))
