#pragma once
#include <shared_mutex>
#include <unordered_set>
#include <cstring>
#include <span>
#include "pbrt/util/hash.h"
#include "pbrt/math/vecmath.h"

namespace pbrt {

    /*
     * Handles storing pointers to unique arrays
     * Based on PBRT
    */
    template <typename T>
    class BufferCache {
        public:
            BufferCache() {}
            ~BufferCache() {
                for (const auto& shardSet: cache) {
                    for (const Buffer& buffer: shardSet) {
                        if (buffer.ptr != nullptr) delete[] buffer.ptr;
                    }
                }
            }
            // Maybe FIXME: PBRT had Allocator element
            const T* LookupOrAdd(std::span<const T> buf) 
            {
                Buffer lookupBuffer(buf.data(), buf.size());

                // Find which shard the lookupBuffer falls in
                int shardIndex = uint32_t(lookupBuffer.hash) >> (32 - logShards);

                // If the buffer exists in the cache, return it
                mutex[shardIndex].lock_shared();
                if (auto iter = cache[shardIndex].find(lookupBuffer); 
                    iter != cache[shardIndex].end()) 
                {
                    const T *ptr = iter->ptr;
                    mutex[shardIndex].unlock_shared();
                    return ptr;
                }
                mutex[shardIndex].unlock_shared();

                // Otherwise, add it to the buffer
                T *ptr = new T[buf.size()];
                std::copy(buf.begin(), buf.end(), ptr);

                // Exclusive lock while inserting into cache
                mutex[shardIndex].lock();

                // Handle the case of another thread adding the buffer first
                //if (auto iter = cache[shardIndex].find(lookupBuffer); iter != cache[shardIndex].end()) {
                    //const T *cachePtr = iter->ptr;
                    //mutex[shardIndex].unlock();
                    //delete[] ptr;
                    //return cachePtr;
                //}

                // If no other thread added it, we add it
                cache[shardIndex].insert(Buffer(ptr, buf.size()));
                mutex[shardIndex].unlock();
                return ptr;
            }
        private:
            struct Buffer {
                const T *ptr = nullptr;
                size_t size = 0, hash;

                Buffer(const T *ptr, size_t size)
                    : ptr(ptr), size(size)
                {
                    hash = HashBuffer(ptr, size);
                }
                bool operator==(const Buffer &b) const {
                    return size == b.size && hash == b.hash &&
                        std::memcmp(ptr, b.ptr, size * sizeof(T)) == 0;
                }
            };

            struct BufferHasher {
                size_t operator()(const Buffer &b) const {
                    return b.hash;
                }
            };

            /* BufferCache allows for multithreading.
             * The cache is borken into 64 shards,
             * each of which has a mutex */
            static constexpr int logShards = 6;
            static constexpr int nShards = 1 << logShards;
            std::shared_mutex mutex[nShards];
            std::unordered_set<Buffer, BufferHasher> cache[nShards];
    };

    void InitBufferCaches();
    extern BufferCache<unsigned int> *uintBufferCache;
    extern BufferCache<Point3f> *point3BufferCache;
}
