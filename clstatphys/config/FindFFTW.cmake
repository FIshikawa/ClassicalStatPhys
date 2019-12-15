if(FFTW_LIBRARIES)
  set (FFTW_FIND_QUIETLY TRUE)
endif(FFTW_LIBRARIES)

SET(_PATHS $ENV{FFTW_ROOT} /opt/local)
find_path(FFTW_INCLUDE_DIR fftw3.h HINTS ${_PATHS} PATH_SUFFIXES include)
find_library(FFTW_LIBRARIES fftw3 HINTS ${_PATHS} PATH_SUFFIXES lib)
UNSET(_PATHS)

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (FFTW DEFAULT_MSG FFTW_LIBRARIES FFTW_INCLUDE_DIR)
mark_as_advanced(FFTW_INCLUDE_DIR, FFTW_LIBRARIES)
