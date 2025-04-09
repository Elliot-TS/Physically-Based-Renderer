add_subdirectory(util)
add_library(pbrt
    shapes.cpp
    aggregates.cpp
    integrator.cpp
    util/parallel.cpp
    util/buffercache.cpp
    math/math.cpp
)

target_include_directories(
    pbrt PUBLIC
    "${PROJECT_SOURCE_DIR}"
)
