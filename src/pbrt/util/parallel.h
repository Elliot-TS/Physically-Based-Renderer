#pragma once
#include <functional>
#include <cstdint>

namespace pbrt{
    void ParallelFor2D 
        (const unsigned int width, 
         const unsigned int height,
         std::function<void(
             const unsigned int startIndex,
             const unsigned int endIndex,
             const uint8_t threadID)>
        );
}
