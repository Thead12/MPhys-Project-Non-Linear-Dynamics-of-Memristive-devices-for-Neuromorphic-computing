# CMAKE file for compiling c++ into python module

cmake_minimum_required(VERSION 3.12)
project(SelectivitySimulator)

# 1 force C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 2 find dependencies
find_package(Python COMPONENTS Interpreter Development REQUIRED)
find_package(pybind11 REQUIRED)
find_package(OpenMP REQUIRED)

# 3 define the module
pybind11_add_module(SelectivityLib SelectivityExperiment.cpp)

# 4 link OpenMP to the module
target_link_libraries(SelectivityLib PRIVATE OpenMP::OpenMP_CXX)

# 5 optimisation flags 
if(MSVC)
    target_compile_options(SelectivityLib PRIVATE /O2 /openmp)
else()
    target_compile_options(SelectivityLib PRIVATE -O3 -fopenmp)
endif()