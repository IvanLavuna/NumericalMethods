cmake_minimum_required(VERSION 3.17)
project(NumericalMethodsLabs)

# variables
set(INCLUDE_PATH /home/ivan/CLionProjects/NumericalMethodsLabs/include)
set(LIB_PATH /home/ivan/CLionProjects/NumericalMethodsLabs/lib)
set(CMAKE_CXX_STANDARD 14)
set(SOURCE_FILES  src/LinearAlgebra.cpp  src/Basic.cpp  src/Matrix.cpp src/main.cpp)
set(TEST_SRC  testing/TestGenerator.cpp testing/IterativeMethods.cpp testing/EigenValuesAndVectors.cpp)
set(LAGRANGE_SRC testing/LagrangePoints.cpp)

include_directories(${INCLUDE_PATH}) # To add include path to cmake(with headers)

link_directories(${LIB_PATH}/lib) # To add library path to project

include_directories(src include testing)

# Project executable
add_executable(NumericalMethodsLabs ${SOURCE_FILES})

# Testing executable
add_executable(Test ${TEST_SRC})

#Creating library
add_library(LinearAlgebra SHARED ${SOURCE_FILES})

# Lagrange executable(Fills points.json file with points)
add_executable(Lagrange ${LAGRANGE_SRC})

target_precompile_headers(NumericalMethodsLabs PRIVATE include/pch.h)

target_link_libraries(Test ${LIB_PATH}/libLinearAlgebra.so)
target_link_libraries(Lagrange  ${LIB_PATH}/libLinearAlgebra.so -ljsoncpp)