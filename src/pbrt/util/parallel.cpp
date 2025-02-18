#include "parallel.h"
#include <vector>
#include <thread>
#include <functional>

void pbrt::ParallelFor2D 
(const unsigned int width, 
 const unsigned int height,
 std::function<void(
     const unsigned int startIndex, 
     const unsigned int endIndex)> forEachTile 
 ){
    const unsigned int threads_per_processor = 1;
    unsigned int processor_count = std::thread::hardware_concurrency();
    processor_count = processor_count == 0 ? 1 : processor_count;

    const unsigned int nThreads = processor_count * threads_per_processor;

    // Tile Size
    const unsigned int totalSpan = width*height;
    const unsigned int tileSpan = totalSpan / nThreads;

    // Create threads
    std::vector<std::thread> threads;
    for (int i = 0; i < nThreads; ++i) {
        const unsigned int span = i == nThreads
            ? totalSpan - tileSpan*i
            : tileSpan;
        const unsigned int startIndex = tileSpan*i;
        threads.push_back(std::thread(
                    forEachTile,
                    startIndex,
                    startIndex + span));
    }

    for (int i = 0; i < threads.size(); ++i) {
        threads[i].join();
    }
}
