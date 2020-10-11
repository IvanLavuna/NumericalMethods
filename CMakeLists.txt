cmake_minimum_required(VERSION 3.17)
project(NumericalMethodsLabs)

set(CMAKE_CXX_STANDARD 14)
set(SOURCE_FILES main.cpp src/SLAE.cpp)

include_directories(src headers)
add_executable(NumericalMethodsLabs ${SOURCE_FILES})



