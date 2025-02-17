cmake_minimum_required(VERSION 3.12)

project(Renderer)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

# Debug or Release
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS "-ltbb")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")


add_subdirectory(pbrt)

find_package(
    OpenCV 4.9.0 REQUIRED
    COMPONENTS core imgproc highgui photo
)
find_package(
    SDL3 REQUIRED CONFIG REQUIRED
    COMPONENTS SDL3-shared
)

add_executable(Renderer WIN32 main.cpp)

include_directories( ${OpenCV_INCLUDE_DIRS} )
target_link_libraries(
    Renderer PUBLIC
    pbrt
    SDL3::SDL3
    ${OpenCV_LIBS}
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
