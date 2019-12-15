set(PROGS
  omp
)
foreach(name ${PROGS})
  set(target_name test_${name})
  add_executable(${target_name} ${name}.cpp)
  set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
  target_link_libraries(${target_name} gtest_main)
  add_test(${target_name} ${name})
endforeach(name)

set(PROGS
  mpi
  hybrid
)
foreach(name ${PROGS})
  set(target_name test_${name})
  add_executable(${target_name} ${name}.cpp)
  set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
  target_link_libraries(${target_name} gtest_main)
  set(MPI_COMMAND mpirun -np 2 ${name})
  add_test(
    NAME ${target_name}
    COMMAND ${MPI_COMMAND}
  )
endforeach(name)
