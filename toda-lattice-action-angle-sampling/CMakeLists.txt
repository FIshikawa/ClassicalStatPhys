cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project(toda-lattice-gge NONE)

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
set(GOOGLETEST_DIR "${PROJECT_SOURCE_DIR}/../clstatphys/tool/googletest/")
add_subdirectory(${GOOGLETEST_DIR} ./build/tool/)
include_directories(${GOOGLETEST_DIR}/include)

# include simulation depended libraries
include_directories(${PROJECT_SOURCE_DIR}/include)

# IO test
include(add_iotest)
enable_testing()

add_subdirectory(src)
add_subdirectory(tests)
