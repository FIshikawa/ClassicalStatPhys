set(PROGS 
  hello_omp
  )
foreach(name ${PROGS})
  add_executable(${name} ${name}.cpp)
  add_iotest(${name})
endforeach(name)

set(PROGS 
  hello_hybrid
  hello_mpi
  )
foreach(name ${PROGS})
  add_executable(${name} ${name}.cpp)
  set(MPI_COMMAND mpirun -np 2 ${name})
  add_test(
    NAME ${name}
    COMMAND ${MPI_COMMAND}
  )
endforeach(name)
