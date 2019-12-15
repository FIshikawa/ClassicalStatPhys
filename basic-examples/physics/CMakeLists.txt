set(PROGS 
  harmonic_oscillator
  todalattice
  fpu
  valence_force_field_model_fixed_end
  exactcalc
  )
foreach(name ${PROGS})
  add_executable(${name} ${name}.cpp)
  target_link_libraries(${name} ${FFTW_LIBRARIES})
  add_iotest(${name})
endforeach(name)
