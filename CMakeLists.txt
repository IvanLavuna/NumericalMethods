cmake_minimum_required(VERSION 3.17)
project(NumericalMethodsLabs)


set(CMAKE_CXX_STANDARD 14)
set(SOURCE_FILES testing/Gauss_Jordan.cpp  src/LinearAlgebra.cpp
        src/Basic.cpp testing/TestGenerator.cpp
        testing/TestGenerator.h src/Matrix.cpp src/Matrix.h)
include_directories(src include testing)

add_executable(NumericalMethodsLabs ${SOURCE_FILES})



