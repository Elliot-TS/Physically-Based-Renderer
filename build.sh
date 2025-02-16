#!/bin/bash

# Set the default build type to Debug
BUILD_TYPE="Debug"

# Check if an argument is provided
if [ $# -gt 0 ]; then
  # Check if the argument is "Debug" or "Release" (case-insensitive)
  case "$1" in
    [Dd]ebug)
      BUILD_TYPE="Debug"
      ;;
    [Rr]elease)
      BUILD_TYPE="Release"
      ;;
    *)  # Invalid argument
      echo "Invalid build type: $1.  Using default: Debug"
      ;;
  esac
fi

cd bin 
cmake ../src -DCMAKE_BUILD_TYPE=$BUILD_TYPE
make
