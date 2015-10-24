cmake_minimum_required(VERSION 2.8.4)

if (POLICY CMP0054)
  cmake_policy(SET CMP0054 NEW)
endif()

if (POLICY CMP0042)
  cmake_policy(SET CMP0042 NEW)
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")
SET(CMAKE_BUILD_TYPE "DEBUG")

PROJECT(mex_functions)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --std=c++11 -fPIC")
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC -DMX_COMPAT_32 -DMATLAB_MEX_FILE")

SET(OCTAVE_FIND_DEBUG True)

find_package(OCTAVE)
if (OCTAVE_FOUND)
  message( STATUS "Found Octave")
	set(MATLAB_LIBRARIES ${OCTAVE_LIBRARIES})
	set(MATLAB_MEX_LIBRARIES ${OCTAVE_OCTINTERP_LIBRARY})
	set(Matlab_MEX_LIBRARY ${OCTAVE_OCTINTERP_LIBRARY})
	set(MATLAB_MX_LIBRARIES ${OCTAVE_OCT_LIB_DIR})
	set(Matlab_MEX_EXTENSION "mex")
	set(MATLAB_MEX_EXTENSION "mex")
else()
  message( STATUS "No Octave! - Looking for MATLAB" )
	find_package( MATLAB COMPONENTS MX_LIBRARY)
	if ( MATLAB_FOUND )
		message( STATUS "Found Matlab at dir:     ${MATLAB_ROOT_DIR}" )
		message( STATUS "Matlab include dir:      ${MATLAB_INCLUDE_DIRS}" )
	endif()
endif()

# Add for tests to be more verbose
#add_definitions(-DDEBUG)

# Add tests for Eigen if found
find_package(EIGEN)
if (EIGEN_FOUND) 
	include_directories(${MATLAB_INCLUDE_DIRS} ${OCTAVE_INCLUDE_DIRS} ${EIGEN_INCLUDE_DIR} )
	set(TESTS test_eig_add test_eig_mult)

	FOREACH(t ${TESTS})
		add_executable(${t} ${t}.cpp)
		target_link_libraries(${t} ${MATLAB_LIBRARIES})
	ENDFOREACH()
	
	set(EIG_MEXES eigen_add eigen_vector_example)

	FOREACH(t ${EIG_MEXES})
		if (OCTAVE_FOUND)
			octave_add_oct(${t} SOURCES ${t}.cpp)
		else()
			add_matlab_mex(NAME ${t} SRC ${t}.cpp	LINK_TO ${MATLAB_LIBRARIES})
		endif()
	ENDFOREACH(t ${EIG_MEXES})
else()
	include_directories(${MATLAB_INCLUDE_DIRS} ${OCTAVE_INCLUDE_DIRS} )
endif(EIGEN_FOUND)


set(MEXES simple_example vector_example)

# This comes from FindMATLAB.cmake
FOREACH(t ${MEXES})
		if (OCTAVE_FOUND)
			octave_add_oct(${t} SOURCES ${t}.cpp)
		else()
			add_matlab_mex(NAME ${t} SRC ${t}.cpp	LINK_TO ${MATLAB_LIBRARIES})
		endif()
ENDFOREACH(t ${MEXES})


# Standalone example with external mex function
add_executable(test_mex test_mex.cpp mex_function.cpp)
target_link_libraries(test_mex ${MATLAB_LIBRARIES})

add_executable(test_vector_add test_vector_add.cpp)
target_link_libraries(test_vector_add ${MATLAB_LIBRARIES})


SET(LIBRARY_OUTPUT_PATH ..)
SET(EXECUTABLE_OUTPUT_PATH ..)
