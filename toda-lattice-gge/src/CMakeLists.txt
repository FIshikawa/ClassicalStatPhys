cmake_minimum_required(VERSION 3.1 FATAL_ERROR)# minimun required version

set(PROGS 
  toda_lattice_periodic_boundary_thermalization_from_normalmode_ensemble
  )
foreach(name ${PROGS})
  add_executable(${name} ${name}.cpp)
  set(test_name "${name}")
  target_link_libraries(${name} ${FFTW_LIBRARIES})
  target_link_libraries(${name} ${ROKKO_LIBRARIES})
  file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/src/test/)
  file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/src/test/${test_name}/)
  set(test_directory "${PROJECT_BINARY_DIR}/src/test/${test_name}/")
  set(MPI_COMMAND mpirun -np 2 ${name})
  add_test(
    NAME ${name}
    COMMAND ${MPI_COMMAND}
  )
endforeach(name)

unset(PROG)
