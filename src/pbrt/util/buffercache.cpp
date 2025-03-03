#include "buffercache.h"

namespace pbrt {
    BufferCache<unsigned int> *uintBufferCache;
BufferCache<Point3f> *point3BufferCache;

void InitBufferCaches() {
    CHECK(uintBufferCache == nullptr);
    uintBufferCache = new BufferCache<unsigned int>;
    point3BufferCache = new BufferCache<Point3f>;
}
}
