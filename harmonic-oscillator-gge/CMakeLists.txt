cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project(harmonic-oscillator-gge NONE)

FILE(GLOB SELFCONFIG RELATIVE ${PROJECT_SOURCE_DIR} config)
if(SELFCONFIG)
  message(STATUS "Find config in " ${PROJECT_SOURCE_DIR})
  list(APPEND CMAKE_MODULE_PATH ${SELFCONFIG})
endif(SELFCONFIG)

set(CLSTATPHYS_DIR "${PROJECT_SOURCE_DIR}/../clstatphys/")
message(STATUS "Install clstatphys from" ${CLSTATPHYS_DIR})
list(APPEND CMAKE_MODULE_PATH ${CLSTATPHYS_DIR}/config)
include(core)
include_directories(${CLSTATPHYS_DIR})

set(GOOGLETEST_VERSION 1)
add_subdirectory(tool/googletest)
include_directories(tool/googletest/include)

# include simulation depended libraries
include_directories(${PROJECT_SOURCE_DIR}/include)

# IO test
include(add_iotest)
enable_testing()

add_subdirectory(src)
add_subdirectory(tests)
