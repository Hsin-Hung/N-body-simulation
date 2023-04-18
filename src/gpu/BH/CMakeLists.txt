cmake_minimum_required(VERSION 3.10)

project(BarnesHut)

find_package(OpenCV REQUIRED)
if(NOT OpenCV_FOUND)
    message(FATAL_ERROR "OpenCV not found")
endif()
include_directories(${OpenCV_INCLUDE_DIRS})

SET(CUDA_SEPARABLE_COMPILATION ON)

# Find CUDA
find_package(CUDA REQUIRED)
if(NOT CUDA_FOUND)
    message(FATAL_ERROR "CUDA not found")
endif()
include_directories(${CUDA_INCLUDE_DIRS})

# set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} ${CUDA_ARCH_FLAGS} -rdc=true")
# Set C++ standard
set(CMAKE_CXX_STANDARD 20)

# Add source files
file(GLOB SOURCES *.cpp *.h *.cu *.cuh)

# Add the executable
cuda_add_executable(BarnesHut ${SOURCES})

# Link GLFW and OpenGL libraries
target_link_libraries(BarnesHut ${OpenCV_LIBS} ${CUDA_TOOLKIT_ROOT_DIR}/lib64/libcudadevrt.a ${CUDA_LIBRARIES})