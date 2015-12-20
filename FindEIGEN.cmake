find_path(EIGEN_INCLUDE_DIR Eigen HINTS "/usr/local/include/eigen3/" "/usr/include/eigen3/" ".." NO_DEFAULT_PATH)
set(EIGEN_INCLUDE_DIRS ${EIGEN_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args (EIGEN DEFAULT_MSG EIGEN_INCLUDE_DIRS)

mark_as_advanced(EIGEN_INCLUDE_DIR)
