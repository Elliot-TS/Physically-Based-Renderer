#pragma once
#include <chrono>

namespace pbrt {
class Benchmark {
 private:
  // using std::chrono::duration;
  // using std::chrono::duration_cast;
  // using std::chrono::high_resolution_clock;
  // using std::chrono::milliseconds;
  std::chrono::time_point<std::chrono::high_resolution_clock>
      start;

 public:
  void Start()
  {
    start = std::chrono::high_resolution_clock::now();
  }
  double GetTime(bool seconds = false)
  {
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> diff =
        end - start;
    return seconds ? diff.count() / 1000.0 : diff.count();
  }
};
}
