set(PROGS 
  normalmode_ensemble
  )
foreach(name ${PROGS})
  add_executable(${name} ${name}.cpp)
  target_link_libraries(${name} ${FFTW_LIBRARIES})
  add_iotest(${name})
endforeach(name)
