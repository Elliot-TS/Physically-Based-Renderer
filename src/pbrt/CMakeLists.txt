add_library(pbrt
    shapes.cpp
    primitive.cpp
)

target_include_directories(
    pbrt PUBLIC
    "${PROJECT_SOURCE_DIR}"
)
