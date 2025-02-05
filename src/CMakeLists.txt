cmake_minimum_required(VERSION 3.12)

project(Renderer)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_subdirectory(pbrt)
#add_library(pbrt
    #pbrt/shapes.cpp
    #pbrt/primitive.cpp
#)

add_executable(Renderer main.cpp)

target_link_libraries(
    Renderer PUBLIC
    pbrt
)

target_include_directories(
    Renderer PUBLIC
    "${PROJECT_BINARY_DIR}"
    "${PROJECT_SOURCE_DIR}"
)

# DEBUG
#add_executable(Debug pbrt/test/sizes.cpp)

#target_link_libraries(
    #Debug PUBLIC
    #pbrt
#)

#target_include_directories(
    #Debug PUBLIC
    #"${PROJECT_BINARY_DIR}"
    ##"${PROJECT_SOURCE_DIR}/pbrt"
#)
