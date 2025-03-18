#include "parallel.h"
#include <functional>
#include <thread>
#include <vector>

void pbrt::ParallelFor2D(
    const unsigned int width,
    const unsigned int height,
    std::function<void(
        const unsigned int startIndex,
        const unsigned int endIndex,
        const uint8_t threadID
    )>
        forEachTile
)
{
  const unsigned int threads_per_processor = 16;
  unsigned int processor_count =
      std::thread::hardware_concurrency();
  processor_count = processor_count == 0 ? 1 : processor_count;

  const unsigned int nThreads =
      processor_count * threads_per_processor;

  // Tile Size
  const unsigned int totalSpan = width * height;
  const unsigned int tileSpan = totalSpan / nThreads;

  // Create threads
  std::vector<std::thread> threads;
  for (int i = 0; i < nThreads; ++i) {
    const unsigned int span =
        i == nThreads ? totalSpan - tileSpan * i : tileSpan;
    const unsigned int startIndex = tileSpan * i;
    threads.push_back(std::thread(
        forEachTile, startIndex, startIndex + span, uint8_t(i)
    ));
  }

  for (int i = 0; i < threads.size(); ++i) {
    threads[i].join();
  }
}
