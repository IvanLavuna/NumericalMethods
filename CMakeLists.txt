cmake_minimum_required(VERSION 3.17)
project(NumericalMethodsLabs)


set(CMAKE_CXX_STANDARD 14)
set(SOURCE_FILES testing/Gauss_Jordan.cpp  src/SLAE.cpp src/LargeDem.cpp
        include/LargeDem.h include/Basic.h src/Basic.cpp)
include_directories(src include testing)

add_executable(NumericalMethodsLabs ${SOURCE_FILES})



