cmake_minimum_required(VERSION 2.8.4)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")
SET(CMAKE_BUILD_TYPE "DEBUG")

PROJECT(mex_functions)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --std=c++11 -fPIC")
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC -DMX_COMPAT_32 -DMATLAB_MEX_FILE")

find_package( MATLAB COMPONENTS MX_LIBRARY)
if ( MATLAB_FOUND )
  message( STATUS "Found Matlab at dir:     ${MATLAB_ROOT_DIR}" )
  message( STATUS "Matlab include dir:      ${MATLAB_INCLUDE_DIRS}" )
else( MATLAB_FOUND )
  message( STATUS "Couldn't find MATLAB - looking for Octave" )
  find_package(OCTAVE)
  if (OCTAVE_FOUND)
    message( STATUS "Found Octave")
		set(MATLAB_INCLUDE_DIRS ${OCTAVE_INCLUDE_DIRS})
		set(MATLAB_LIBRARIES ${OCTAVE_LIBRARIES})
		set(MATLAB_MEX_LIBRARIES ${OCTAVE_OCTINTERP_LIBRARY})
		set(MATLAB_MX_LIBRARIES ${OCTAVE_OCT_LIB_DIR})
		set(MATLAB_MEX_EXTENSION "mex")
  else()
    message( STATUS "No Octave!")
  endif( OCTAVE_FOUND )
endif(MATLAB_FOUND)

find_package(EIGEN)

include_directories(${MATLAB_INCLUDE_DIRS} ${EIGEN_INCLUDE_DIRS} "../.." "." )

if (POLICY CMP0054)
  cmake_policy(SET CMP0054 NEW)
endif()

set(MEXES example_mex_function)

FOREACH(t ${MEXES})
  add_matlab_mex(NAME ${t} 
		SRC ${t}.cxx 
		LINK_TO ${MATLAB_LIBRARIES})
ENDFOREACH(t ${MEXES})

add_executable(test_mex test_mex.cpp)
target_link_libraries(test_mex ${MATLAB_LIBRARIES})

SET(LIBRARY_OUTPUT_PATH ..)
SET(EXECUTABLE_OUTPUT_PATH ..)
