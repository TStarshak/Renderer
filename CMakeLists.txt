# CMakeList.txt : CMake project for OctTreeTracer, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

project ("OctTreeTracer")

# Add source to this project's executable.
add_executable (OctTreeTracer "main.cpp"  "modeldata.h" "vec3.cpp" "vec3.h" "matrix44.cpp" "matrix44.h"  "boundingbox.cpp" "boundingbox.h"   "mesh.cpp" "mesh.h")

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET OctTreeTracer PROPERTY CXX_STANDARD 20)
endif()

# TODO: Add tests and install targets if needed.
