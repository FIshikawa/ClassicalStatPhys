set(PROGS 
  1spin_rk4
  1spin_velver
  )
foreach(name ${PROGS})
  add_executable(${name} ${name}.cpp)
  add_iotest(${name})
endforeach(name)
