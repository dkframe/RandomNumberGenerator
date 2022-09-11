# CMakeList.txt : CMake project for RandomNumberGenerator, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

project ("RandomNumberGenerator")

# Add source to this project's executable.
add_executable (RandomNumberGenerator "RandomNumberGenerator.cpp" "RandomNumberGenerator.h")

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET RandomNumberGenerator PROPERTY CXX_STANDARD 20)
endif()

# TODO: Add tests and install targets if needed.
