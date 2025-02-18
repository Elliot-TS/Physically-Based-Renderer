add_subdirectory(util)
add_library(pbrt
    shapes.cpp
    primitive.cpp
    integrator.cpp
    util/parallel.cpp
)

target_include_directories(
    pbrt PUBLIC
    "${PROJECT_SOURCE_DIR}"
)
