cmake_minimum_required(VERSION 3.10)

project(BarnesHut)

find_package(OpenCV REQUIRED)
if(NOT OpenCV_FOUND)
    message(FATAL_ERROR "OpenCV not found")
endif()
include_directories(${OpenCV_INCLUDE_DIRS})

# Find CUDA
find_package(CUDA REQUIRED)
if(NOT CUDA_FOUND)
    message(FATAL_ERROR "CUDA not found")
endif()
include_directories(${CUDA_INCLUDE_DIRS})

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)

# Add source files
set(SOURCES 
    main.cpp
    constants.h
    barnesHut_kernel.cu 
    barnesHut_kernel.cuh 
    barnesHutCuda.cu
    barnesHutCuda.cuh 
)

# Add the executable
cuda_add_executable(BarnesHut ${SOURCES})

# Link GLFW and OpenGL libraries
target_link_libraries(BarnesHut ${OpenCV_LIBS})