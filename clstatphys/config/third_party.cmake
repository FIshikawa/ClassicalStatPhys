# find FFTW library
include_directories($ENV{FFTW_ROOT}/include)
find_package(FFTW REQUIRED)

# find rokko library
find_package(Rokko REQUIRED HINTS ${ROKKO_ROOT_DIR} $ENV{ROKKO_ROOT})
message(STATUS "Found Rokko: ${ROKKO_ROOT_DIR}")
include(${ROKKO_USE_FILE})

#find OpenMP
find_package(OpenMP)
if (OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

#find MPI
find_package(MPI)
if (MPI_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MPI_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_EXE_LINKER_FLAGS}")
endif()

