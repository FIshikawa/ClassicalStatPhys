cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project(clstatphys NONE)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/config)
include(core)
include_directories(${PROJECT_SOURCE_DIR})

set(GOOGLETEST_VERSION 1)
add_subdirectory(tool/googletest)
include_directories(tool/googletest/include)

enable_testing()

add_subdirectory(tests)
