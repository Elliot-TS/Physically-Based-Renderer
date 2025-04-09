add_subdirectory(util)
add_library(pbrt
    shapes.cpp
    aggregates.cpp
    integrator.cpp
    math/math.cpp
    util/buffercache.cpp
    util/parallel.cpp
    util/transform.cpp
)

target_include_directories(
    pbrt PUBLIC
    "${PROJECT_SOURCE_DIR}"
)
