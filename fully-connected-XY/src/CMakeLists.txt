cmake_minimum_required(VERSION 3.1 FATAL_ERROR)# minimun required version

set(PROGS 
  clXYmodel
  clXYmodelCriticalMD 
  clXYmodelCritical 
  clXYmodelNonEq
  clXYmodelNonEqwithBath
  )
set(LATTICE 
  _Square 
  _Fully 
  _Cube 
  _Tesseract
  )
foreach(name ${PROGS}) #loop name in PROGS
  foreach(lattice ${LATTICE})
    set(target_name "${name}${lattice}")
    set(test_name "${target_name}")
    add_executable(${target_name} ${target_name}.cpp)
    target_link_libraries(${target_name} ${FFTW_LIBRARIES})#link include header file to exe file.
    target_link_libraries(${target_name} ${ALPSCore_LIBRARIES})
    file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/src/test/)
    file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/src/test/${test_name}/)
    set(test_directory "${PROJECT_BINARY_DIR}/src/test/${test_name}/")
    add_test(
      NAME ${test_name} 
      COMMAND ${target_name} 
      WORKING_DIRECTORY ${test_directory}
    ) 
  endforeach(lattice)
endforeach(name)
unset(PROG)

set(PROGS 
  #clXYmodel_FullyParallelized 
  #clXYmodelWangLandau_Fully
  )
foreach(name ${PROGS}) #loop name in PROGS
  add_executable(${name} ${name}.cpp)
  set(test_name "${name}")
  target_link_libraries(${name} ${FFTW_LIBRARIES})#link include header file to exe file.
  target_link_libraries(${name} ${ALPSCore_LIBRARIES})
  file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/src/test/)
  file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/src/test/${test_name}/)
  set(test_directory "${PROJECT_BINARY_DIR}/src/test/${test_name}/")
  add_test(
    NAME ${test_name} 
    COMMAND ${name} 
    WORKING_DIRECTORY ${test_directory}$
  ) 

endforeach(name)

unset(PROG)
