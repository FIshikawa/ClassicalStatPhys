cmake_minimum_required(VERSION 3.1 FATAL_ERROR)# minimun required version

FILE(GLOB PROGS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)
foreach(progs ${PROGS})
  string(REGEX REPLACE ".cpp" "" name ${progs}) 
  add_executable(${name} ${name}.cpp)
  set(test_name "${name}")
  target_link_libraries(${name} ${FFTW_LIBRARIES})
  target_link_libraries(${name} ${ROKKO_LIBRARIES})
  set(MPI_COMMAND mpirun -np 1 ${name})
  add_test(
    NAME ${name}
    COMMAND ${MPI_COMMAND}
  )
endforeach(progs)

unset(PROG)
