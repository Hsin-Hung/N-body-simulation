cmake_minimum_required(VERSION 3.21)

# set the project name and version
project(NBody)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

message(STATUS "system: ${CMAKE_SYSTEM}")
message(STATUS "build_type: ${CMAKE_BUILD_TYPE}")
message(STATUS "paths:")
message(STATUS " - ${PROJECT_SOURCE_DIR}")
message(STATUS " - ${CMAKE_INSTALL_PREFIX}")
message(STATUS " - ${CMAKE_SOURCE_DIR}")
message(STATUS " - ${CMAKE_CURRENT_SOURCE_DIR}")
message(STATUS)
message(STATUS "compiler:")
message(STATUS " - ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS " - ${CMAKE_CXX_COMPILER}")
message(STATUS " - ${CMAKE_CXX_FLAGS}")
message(STATUS)
message(STATUS "cmake:")
message(STATUS " - ${CMAKE_VERSION}")
message(STATUS " - ${CMAKE_COMMAND}")
message(STATUS " - ${CMAKE_TOOLCHAIN_FILE}")
message(STATUS " - ${CMAKE_GENERATOR}")
message(STATUS " - ${CMAKE_BUILD_TOOL}")
message(STATUS)

find_package(OpenGL REQUIRED)

if(NOT OPENGL_FOUND)
    message(FATAL_ERROR "OpenGL not found")
endif()

find_package(glfw3 3.3 REQUIRED)

if(NOT glfw3_FOUND)
    message(FATAL_ERROR "GLFW not found")
endif()

file(GLOB SRC_FILES ${PROJECT_SOURCE_DIR}/src/cpu/*.cpp)

add_executable(NBody ${PROJECT_SOURCE_DIR}/src/main.cpp ${SRC_FILES})
target_include_directories(NBody PRIVATE ${PROJECT_SOURCE_DIR}/src/cpu ${GL_INCLUDE_DIRS} ${GLFW_INCLUDE_DIRS})
target_link_libraries(NBody OpenGL::GL glfw)