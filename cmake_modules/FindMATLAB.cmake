#.rst:
# FindMatlab
# ----------------
#
# Finds Matlab installations and provides Matlab tools and libraries to cmake.
# 
# This package first intention is finding the libraries associated with Matlab in order
# to be able to compile Matlab extensions (mex files). It can also be used to run unit test on these mex extensions,
# and run Matlab.
#
# The variable ``Matlab_USER_ROOT`` may be specified in order to give the path of the desired Matlab version.  
# Otherwise, the behaviour is platform dependant:
#
# - on Windows, the installed versions of Matlab are retrieved from the Windows registry
# - on Mac, the installed versions of Matlab are given by the Matlab pathes in /Application
# - on Unix, the desired Matlab should be accessible from the PATH.
#
# Additional information is provided when ``Matlab_FIND_DEBUG`` is set.
# When a Matlab binary is found automatically and the ``Matlab_VERSION`` is not given, then the
# version is queried from Matlab directly. On Windows, it can make a window running Matlab appear.
#
# The mapping of the release names and the version of Matlab is performed by defining pairs (name, version). 
# The variable ``Matlab_ADDITIONAL_VERSIONS`` may be provided in order to handle additional versions, in the following form:
#
# ::
#
#    set(Matlab_ADDITIONAL_VERSIONS 
#       "release_name1" "corresponding_version1"
#       "release_name2" "corresponding_version2"
#       ...
#       )
# 
# such as 
# ::
# 
#    set(Matlab_ADDITIONAL_VERSIONS
#        "R2013b" "8.2"
#        "R2013a" "8.1"
#        "R2012b" "8.0")
#
# Starting Matlab R2013b+, it is possible to use Matlab unit test framework. The syntax is the following
# ::
# 
#    add_matlab_unit_test(
#       NAME <name>
#       UNITTEST_FILE matlab_file_containing_unittest.m
#       [TIMEOUT timeout]
#       [ADDITIONAL_PATH path1 path2 ...]
#       [TEST_ARGS additionnal_arguments1 ...]
#    )
#  
# in particular, the arguments accepted by add_test can be passed through TEST_ARGS (eg. CONFIGURATION <config> ...).
#
# Defined variables
# -----------------
# * ``Matlab_FOUND`` true if the Matlab installation is found.
# * ``Matlab_USER_ROOT`` the root of the Matlab installation. This is given by the user.
# * ``Matlab_ROOT_DIR`` the root of the Matlab installation determined by the FindMatlab module.
# * ``Matlab_VERSION_STRING`` the version of the Matlab installation
# * ``Matlab_PROGRAM`` the Matlab binary program. Available only if the component ``MAIN_PROGRAM`` is asked
# * ``Matlab_INCLUDE_DIRS`` the path of the Matlab libraries headers
# * ``Matlab_MEX_LIBRARY`` library for mex
# * ``Matlab_MX_LIBRARY`` mx library of Matlab (arrays). Available only if the component ``MX_LIBRARY`` is asked
# * ``Matlab_ENG_LIBRARY`` Matlab engine library. Available only if the component ``ENG_LIBRARY`` is asked
# * ``Matlab_LIBRARIES`` the whole set of libraries of Matlab
# * ``Matlab_MEX_COMPILER`` the mex compiler of Matlab. Currently not used internally. Available only if the component ``MEX_COMPILER`` is asked
# * ``Matlab_MEX_EXTENSION`` the extension of the mex files for the current platform (given by Matlab).
#
#
# Defined macros
# --------------
# * ``matlab_get_version_from_release_name`` returns the version from the release name
# * ``matlab_get_release_name_from_version`` returns the release name from the Matlab version
#
# Defined functions
# --------------
# * ``add_matlab_mex`` adds a target compiling a MEX file.
# * ``add_matlab_unit_test`` adds a Matlab unit test file as a test to the project.
# * ``matlab_extract_all_installed_versions_from_registry`` parses the registry for all Matlab versions. Available on Windows only. 
#   The part of the registry parsed is dependent on the host processor 
# * ``matlab_get_all_valid_matlab_roots_from_registry`` returns all the possible Matlab paths, according to a previously given list. Only the
#   existing/accessible paths are kept. This is mainly useful for the "brute force" search of Matlab installation.
# * ``matlab_get_mex_suffix`` returns the suffix to be used for the mex files (platform/architecture dependant)
# * ``matlab_get_version_from_matlab_run`` returns the version of Matlab, given the full directory of the Matlab program.
# 
# Future work
# -----------
# - win32:an additional variable telling that the registry is x86 or x64, maybe depending on the target build.

#=============================================================================
# Copyright 2009-2014 Raffi Enficiaud
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distribute this file outside of CMake, substitute the full
#  License text for the above reference.)

include(CheckCXXCompilerFlag)

# The currently supported versions. Other version can be added by the user by providing Matlab_ADDITIONAL_VERSIONS
if(NOT Matlab_ADDITIONAL_VERSIONS)
  set(Matlab_ADDITIONAL_VERSIONS)
endif()

set(Matlab_FIND_DEBUG FALSE)
set(Matlab_VERSIONS_MAPPING
        "R2014b" "8.4"
  "R2013b" "8.2"
  "R2013a" "8.1"
  "R2012b" "8.0"
  "R2012a" "7.14"

  "R2011b" "7.13"
  "R2011a" "7.12"
  "R2010b" "7.11"
  
  ${Matlab_ADDITIONAL_VERSIONS}
  )


# temporary folder for all Matlab runs
set(_matlab_temporary_folder ${CMAKE_BINARY_DIR}/Matlab)
  
if(NOT EXISTS ${_matlab_temporary_folder})
  file(MAKE_DIRECTORY ${_matlab_temporary_folder})
endif()

# get the version of Matlab (17.58) from a release name (R2017k)
macro (matlab_get_version_from_release_name release_name version_name)
  list(FIND Matlab_VERSIONS_MAPPING ${release_name} _index)
  if(${_index} EQUAL -1)
    message(WARNING "The release name ${release_name} is not registered")
  endif()
  math(EXPR _index "${_index}+1")
  list(GET Matlab_VERSIONS_MAPPING ${_index} _version)
  set(${version_name} ${_version})
  unset(_index)
  unset(_version)
endmacro(matlab_get_version_from_release_name)


# get the release name (R2017k) from the version of Matlab (17.58)
macro (matlab_get_release_name_from_version version release_name)
  list(FIND Matlab_VERSIONS_MAPPING ${version} _index)
  if(${_index} EQUAL -1)
    message(WARNING "The version ${version} is not registered")
  endif()
  math(EXPR _index "${_index}-1")
  list(GET Matlab_VERSIONS_MAPPING ${_index} _release)
  set(${release_name} ${_release})
  unset(_release)
  unset(_index)
endmacro(matlab_get_release_name_from_version)


# extracts all the supported release names (R2017k...) of Matlab
macro(matlab_get_supported_releases list_releases)
  list(LENGTH Matlab_VERSIONS_MAPPING versions_length)
  math(EXPR versions_length "${versions_length}-1")
  set(${list_releases})
  foreach(matlab_release RANGE 0 ${versions_length} 2)
    list(GET Matlab_VERSIONS_MAPPING ${matlab_release} current)
    list(APPEND ${list_releases} ${current})
  endforeach(matlab_release)
endmacro(matlab_get_supported_releases)

# extracts all the supported versions of Matlab
macro(matlab_get_supported_versions list_versions)
  list(LENGTH Matlab_VERSIONS_MAPPING versions_length)
  set(${list_versions})
  foreach(matlab_version RANGE 1 ${versions_length} 2)
    list(GET Matlab_VERSIONS_MAPPING ${matlab_version} current)
    list(APPEND ${list_versions} ${current})
  endforeach(matlab_version)
endmacro(matlab_get_supported_versions)



# This function parses the registry and founds the Matlab versions that are "really" installed
# the other approach is to use a brute force search
# set win64 to TRUE if the 64 bit version of Matlab should be looked for
# the returned list contains all versions under HKLM\\SOFTWARE\\Mathworks\\Matlab or an empty list in case an error occurred (or nothing found)
#
# Only the version is provided: no path, no test for existence
function(matlab_extract_all_installed_versions_from_registry win64 matlab_versions)
  
  if(NOT CMAKE_HOST_WIN32)
    message(FATAL_ERROR "This macro can only be called by a windows host (call to reg.exe")
  endif()
  
  #message(STATUS "System processor ${CMAKE_HOST_SYSTEM_PROCESSOR}")
  
  
  
  # list the keys under HKEY_LOCAL_MACHINE\SOFTWARE\mathworks but the call to reg does not work
  # from cmake, curiously, as is. The command provides the desired result under the command line though.
  # Fix: this is because "/reg:64" should appended to the command, otherwise it gets on the 32 bits software key (curiously again)
  find_program(Matlab_REG_EXE_LOCATION "reg")
  file(TO_NATIVE_PATH ${Matlab_REG_EXE_LOCATION} Matlab_REG_EXE_LOCATION)
  
  # if reg.exe is not found, then it is impossible to use this method.
  if(NOT Matlab_REG_EXE_LOCATION)
    if(Matlab_FIND_DEBUG)
      message(WARNING "[Matlab] reg.exe not found")
    endif()
    set(${matlab_versions} "" PARENT_SCOPE)
    return()
  endif()
  
  
  
  if(${win64} AND ${CMAKE_HOST_SYSTEM_PROCESSOR} MATCHES "64")
    set(APPEND_REG "/reg:64")
  else()
    set(APPEND_REG "/reg:32")
  endif()
  
  # /reg:64 should be added on 64 bits capable OSs in order to enable the redirection of 64 bits applications
  execute_process(
    COMMAND ${Matlab_REG_EXE_LOCATION} query HKEY_LOCAL_MACHINE\\SOFTWARE\\Mathworks\\Matlab /f * /k ${APPEND_REG}
    RESULT_VARIABLE resultMatlab
    OUTPUT_VARIABLE varMatlab)
  #message("Matlabs = ${varMatlab} | ${resultMatlab}")
  
  
  set(matlabs_from_registry)
  if(${resultMatlab} EQUAL 0)
    # just for tests
    #set(varMatlab "${varMatlab} HKEY_LOCAL_MACHINE\\SOFTWARE\\Mathworks\\Matlab\\7.47 HKEY_LOCAL_MACHINE\\SOFTWARE\\Mathworks\\Matlab\\9")
    string(
      REGEX MATCHALL "Matlab\\\\([0-9]+(\\.[0-9]+)?)"
      matlab_versions_regex ${varMatlab})
    #message(STATUS "regex = ${matlab_versions_regex}")
    foreach(match IN LISTS matlab_versions_regex)
      string(
        REGEX MATCH "Matlab\\\\(([0-9]+)(\\.([0-9]+))?)"
        current_match ${match})
      #message(STATUS "current match is ${CMAKE_MATCH_0} ${CMAKE_MATCH_1} ${CMAKE_MATCH_2} ${CMAKE_MATCH_3} ${CMAKE_MATCH_4}")
      
      set(_matlab_current_version ${CMAKE_MATCH_1})
      set(current_matlab_version_major ${CMAKE_MATCH_2})
      set(current_matlab_version_minor ${CMAKE_MATCH_4})
      if(NOT current_matlab_version_minor)
        set(current_matlab_version_minor "0")
      endif()
      
      #message(STATUS "Matlab in registry ${current_matlab_version_major}.${current_matlab_version_minor}")

      list(APPEND matlabs_from_registry ${_matlab_current_version})
      unset(_matlab_current_version)
    endforeach(match)
    
  endif()
  
  set(${matlab_versions} ${matlabs_from_registry} PARENT_SCOPE)

endfunction(matlab_extract_all_installed_versions_from_registry)

macro(extract_matlab_versions_from_registry_brute_force matlab_versions)
  # get the supported versions
  set(matlab_supported_versions)
  matlab_get_supported_versions(matlab_supported_versions)
  
  
  # this is a manual population of the versions we want to look for
  # this can be done as is, but preferably with the call to 
  # matlab_get_supported_versions and variable 
  
  # populating the versions we want to look for
  # set(matlab_supported_versions)
  
  # # Matlab 7
  # set(matlab_major 7)
  # foreach(current_matlab_minor RANGE 4 20)
    # list(APPEND matlab_supported_versions "${matlab_major}.${current_matlab_minor}")
  # endforeach(current_matlab_minor)

  # # Matlab 8
  # set(matlab_major 8)
  # foreach(current_matlab_minor RANGE 0 5)
    # list(APPEND matlab_supported_versions "${matlab_major}.${current_matlab_minor}")
  # endforeach(current_matlab_minor)
  
  # # taking into account the possible additional versions provided by the user
  # if(DEFINED Matlab_ADDITIONAL_VERSIONS)
    # list(APPEND matlab_supported_versions Matlab_ADDITIONAL_VERSIONS)
  # endif()
  
  
  # we order from more recent to older
  list(REMOVE_DUPLICATES matlab_supported_versions)
  list(SORT matlab_supported_versions)
  list(REVERSE matlab_supported_versions)
  
  
  set(${matlab_versions} ${matlab_supported_versions})
  

endmacro(extract_matlab_versions_from_registry_brute_force)


# populates the Matlab root with valid versions of Matlab. The matlab_versions comes either from
# extract_matlab_versions_from_registry_brute_force or matlab_extract_all_installed_versions_from_registry.
# The returned matlab_roots is organized in pairs version_number,matlab_root_path.  
function(matlab_get_all_valid_matlab_roots_from_registry matlab_versions matlab_roots)
  
  set(_matlab_roots_list )
  foreach(_matlab_current_version ${matlab_versions})
    get_filename_component(
      current_Matlab_ROOT
      "[HKEY_LOCAL_MACHINE\\SOFTWARE\\MathWorks\\Matlab\\${_matlab_current_version};MatlabROOT]"
      ABSOLUTE)
      
    if(EXISTS ${current_Matlab_ROOT})
      list(APPEND _matlab_roots_list ${_matlab_current_version} ${current_Matlab_ROOT})
    endif()
  
  endforeach(_matlab_current_version)
  unset(_matlab_current_version)
  set(${matlab_roots} ${_matlab_roots_list} PARENT_SCOPE)
  unset(_matlab_roots_list)
endfunction(matlab_get_all_valid_matlab_roots_from_registry)


# returns the extension of the mex files (the suffixes).
# This function should not be called before the appropriate Matlab root has been found.
function(matlab_get_mex_suffix matlab_root mex_suffix)

  # todo setup the extension properly. Currently I do not know if this is sufficient for all win32 distributions.
  # there is also CMAKE_EXECUTABLE_SUFFIX that could be tweaked
  set(mexext_suffix "")
  if(WIN32)
    list(APPEND mexext_suffix ".bat")
  endif()

  # we first try without suffix, since cmake does not understand a list with one empty string element
  find_program(
    Matlab_MEXEXTENSIONS_PROG
    "mexext"
    PATHS ${matlab_root}/bin
    DOC "Matlab MEX extension provider"
    NO_DEFAULT_PATH
  )

  foreach(current_mexext_suffix IN LISTS mexext_suffix)
    if(NOT DEFINED Matlab_MEXEXTENSIONS_PROG OR NOT Matlab_MEXEXTENSIONS_PROG)
      # this call should populate the cache automatically
      find_program(
        Matlab_MEXEXTENSIONS_PROG
        "mexext${current_mexext_suffix}"
        PATHS ${matlab_root}/bin
        DOC "Matlab MEX extension provider"
        NO_DEFAULT_PATH
      )
    endif()
  endforeach(current_mexext_suffix)
  
  
  # the program has been found?
  if((NOT Matlab_MEXEXTENSIONS_PROG) OR (NOT EXISTS ${Matlab_MEXEXTENSIONS_PROG}))
    if(Matlab_FIND_DEBUG)
      message(WARNING "[Matlab] Cannot found mexext program. Matlab root is ${matlab_root}")
    endif()
    return()
  endif()

  set(_matlab_mex_extension)

  
  execute_process(
    COMMAND ${Matlab_MEXEXTENSIONS_PROG} 
    OUTPUT_VARIABLE _matlab_mex_extension)
  string(STRIP ${_matlab_mex_extension} _matlab_mex_extension)

  set(${mex_suffix} ${_matlab_mex_extension} PARENT_SCOPE)
endfunction(matlab_get_mex_suffix)



# .. command:: matlab_get_version_from_matlab_run
#
#   This function runs Matlab specified on arguments and extracts its version.
#   
#   matlab_get_version_from_matlab_run(matlab_binary_path matlab_version)
function(matlab_get_version_from_matlab_run matlab_binary_program matlab_list_versions)

  set(${matlab_list_versions} "" PARENT_SCOPE)

  
  if(Matlab_FIND_DEBUG)
    message(STATUS "[Matlab] Determining the version of Matlab from ${matlab_binary_program}")
  endif()

  if(EXISTS ${_matlab_temporary_folder}/matlabVersionLog.cmaketmp)
    if(Matlab_FIND_DEBUG)
      message(STATUS "[Matlab] Removing previous ${_matlab_temporary_folder}/matlabVersionLog.cmaketmp file")
    endif()
    file(REMOVE ${_matlab_temporary_folder}/matlabVersionLog.cmaketmp)
  endif()

  
  # the log file is needed since on windows the command executes in a new window and it is not possible 
  # to get back the answer of Matlab
  # the -wait command is needed on windows, otherwise the call returns immediately after the program launches itself.
  if(WIN32)
    set(_matlab_additional_commands "-wait")
  endif()
  
  # timeout set to 30 seconds, in case it does not start
  # note as said before OUTPUT_VARIABLE cannot be used in a platform independent manner
  # however, not setting it would flush the output of Matlab in the current console (unix variant)
  execute_process(
    COMMAND ${matlab_binary_program} -nosplash -nojvm ${_matlab_additional_commands} -logfile ${_matlab_temporary_folder}/matlabVersionLog.cmaketmp -nodesktop -nodisplay -r "version, exit" 
    OUTPUT_VARIABLE _matlab_version_from_cmd_dummy
    RESULT_VARIABLE _matlab_result_version_call
    TIMEOUT 30
    )
  
  
  if(${_matlab_result_version_call})
    if(Matlab_FIND_DEBUG)
      message(WARNING "[Matlab] Unable to determine the version of Matlab. Matlab call returned with error ${_matlab_result_version_call}.")
    endif()
    return()
  elseif(NOT EXISTS ${_matlab_temporary_folder}/matlabVersionLog.cmaketmp)
    if(Matlab_FIND_DEBUG)
      message(WARNING "[Matlab] Unable to determine the version of Matlab. The log file does not exist.")
    endif()
    return()
  endif()

  # if successful, read back the log
  file(READ ${_matlab_temporary_folder}/matlabVersionLog.cmaketmp _matlab_version_from_cmd)
  file(REMOVE ${_matlab_temporary_folder}/matlabVersionLog.cmaketmp)

  set(index -1)
  string(FIND ${_matlab_version_from_cmd} "ans" index)
  if(index EQUAL -1)
    
    if(Matlab_FIND_DEBUG)
      message(WARNING "[Matlab] Cannot find the version of Matlab returned by the run.")
    endif()
    
  else()
    set(matlab_list_of_all_versions_tmp)
    
    string(SUBSTRING ${_matlab_version_from_cmd} ${index} -1 substring_ans)
    string(
      REGEX MATCHALL "ans[\r\n\t ]*=[\r\n\t ]*([0-9]+(\\.[0-9]+)?)"
      matlab_versions_regex 
      ${substring_ans})
    foreach(match IN LISTS matlab_versions_regex)
      string(
        REGEX MATCH "ans[\r\n\t ]*=[\r\n\t ]*(([0-9]+)(\\.([0-9]+))?)"
        current_match ${match})
      
      list(APPEND matlab_list_of_all_versions_tmp ${CMAKE_MATCH_1})
    endforeach()
    list(REMOVE_DUPLICATES matlab_list_of_all_versions_tmp)
    set(${matlab_list_versions} ${matlab_list_of_all_versions_tmp} PARENT_SCOPE)
    
  endif()
    
endfunction(matlab_get_version_from_matlab_run)


# .. command:: add_matlab_unit_test
#
#   Adds a Matlab unit test. 
#   The unit test uses the Matlab unittest framework, which is available starting Matlab 2013b+. The function expects
#   one Matlab unit test file to be given.
#   The function arguments are:
#   * `NAME` name of the unittest in ctest.
#   * `UNITTEST_FILE` the matlab unittest file. Its path will be automatically added to the Matlab path.
#   * `TIMEOUT` the test timeout in seconds (optional, defaults to 180 seconds).
#   * `ADDITIONAL_PATH` a list of paths to add to the Matlab path prior to running the unit test.
#   * `Matlab_ADDITIONAL_STARTUP_OPTIONS` a list of additional option in order to run Matlab from the command line. 
#   * `TEST_ARGS` Additional options provided to the add_test command. These options are added to the default 
#      options (eg. "CONFIGURATIONS Release")
#
function(add_matlab_unit_test)
  set(options )
  set(oneValueArgs NAME UNITTEST_FILE TIMEOUT)
  set(multiValueArgs ADDITIONAL_PATH Matlab_ADDITIONAL_STARTUP_OPTIONS TEST_ARGS)
  
  set(prefix _matlab_unittest_prefix)
  cmake_parse_arguments(${prefix} "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )
  
  if(NOT ${prefix}_NAME)
    message(FATAL_ERROR "[Matlab] The Matlab test name cannot be empty")
  endif()

  add_test(NAME ${${prefix}_NAME}
           COMMAND ${CMAKE_COMMAND} 
            -Dtest_name=${${prefix}_NAME}
            -Dadditional_paths=${${prefix}_ADDITIONAL_PATH}
            -Dtest_timeout=${${prefix}_TIMEOUT}
            -Doutput_directory=${_matlab_temporary_folder}
            -DMatlab_PROGRAM=${Matlab_PROGRAM}
            -DMatlab_ADDITIONNAL_STARTUP_OPTIONS=${${prefix}_Matlab_ADDITIONAL_STARTUP_OPTIONS}
            -Dunittest_file_to_run=${${prefix}_UNITTEST_FILE}
            -P ${CMAKE_ROOT}/Modules/FindMatlab_TestsRedirect.cmake   ## TODO not sure about the path, to be tested
           ${${prefix}_TEST_ARGS}
           ${${prefix}_UNPARSED_ARGUMENTS}
           )
endfunction(add_matlab_unit_test)


# .. command:: add_matlab_mex
#
#   Adds a Matlab MEX file. 
#   ::
#       add_matlab_mex(
#          NAME <name>
#          SRC src1 [src2 ...]
#          [REDUCE_VISIBILITY]
#          [OUTPUT_NAME output_name]
#          [DOCUMENTATION documentation.m]
#          [LINK_TO target1 target2 ...]
#       )
#
#   `NAME` name of the target. `SRC` is a list of source files. `LINK_TO` is a list of additional link dependencies. 
#   The target already links to libmex and libmx.
#   If `OUTPUT_NAME` is given, it overrides the default name defined as being the name of the target without any prefix and 
#   with `${Matlab_MEX_EXTENSION}` suffix.
#   Use the option `REDUCE_VISIBILITY` to hide every symbols inside the compiled MEX file. This is useful in case their is a 
#   symbol collision between the libraries shipped with Matlab, and the libraries to which the target is linking with. This 
#   option is ignored on DLL platforms.
#   If `DOCUMENTATION` is given, the file documentation.m is considered as being the documentation file for the MEX file. This
#   file is copied into the same folder with the same name as the final mex file, and with extension .m. In that case, typing
#   `help <name>` should print the documentation contained in this file.
#
function(add_matlab_mex )

  set(options REDUCE_VISIBILITY)
  set(oneValueArgs NAME DOCUMENTATION OUTPUT_NAME)
  set(multiValueArgs LINK_TO SRC)
  
  set(prefix _matlab_addmex_prefix)
  cmake_parse_arguments(${prefix} "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )
  
    if (Matlab_FIND_DEBUG)
        message(STATUS "[Matlab] ${oneValueArgs} && ${multiValueArgs}" && "${prefix}_NAME" && ${prefix}_OUTPUT_NAME)
    endif ()
  if(NOT ${prefix}_NAME)
    message(FATAL_ERROR "[Matlab] The MEX file name cannot be empty")
  endif()
  
  if(NOT ${prefix}_OUTPUT_NAME)
    set(${prefix}_OUTPUT_NAME ${${prefix}_NAME})
  endif()
  
    if (Matlab_FIND_DEBUG)
        message(STATUS "[Matlab] output_name ${${prefix}_NAME}")
    endif ()
  add_library(${${prefix}_NAME}
    SHARED 
      ${${prefix}_SRC}
      ${${prefix}_DOCUMENTATION}
      ${${prefix}_UNPARSED_ARGUMENTS})
    if (Matlab_FIND_DEBUG)
        message(STATUS "[Matlab] output_name ${${prefix}_SRC} &&
            ${${prefix}_DOCUMENTATION} &&
            ${${prefix}_UNPARSED_ARGUMENTS}")
    endif ()
  target_include_directories(${${prefix}_NAME} PRIVATE ${Matlab_INCLUDE_DIRS})
  
    if (Matlab_FIND_DEBUG)
        message(STATUS "libraries ${Matlab_MEX_LIBRARY} ${Matlab_MX_LIBRARY} ${${prefix}_LINK_TO}")
    endif ()
  target_link_libraries(${${prefix}_NAME} ${Matlab_MEX_LIBRARY} ${Matlab_MX_LIBRARY} ${${prefix}_LINK_TO})
  set_target_properties(${${prefix}_NAME}
      PROPERTIES 
        PREFIX ""
        OUTPUT_NAME ${${prefix}_OUTPUT_NAME}
        SUFFIX ".${Matlab_MEX_EXTENSION}")


  if(NOT ${${prefix}_DOCUMENTATION} STREQUAL "")
    get_target_property(output_name ${${prefix}_NAME} OUTPUT_NAME)
    add_custom_command(
      TARGET ${${prefix}_NAME}
      PRE_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different ${${prefix}_DOCUMENTATION} $<TARGET_FILE_DIR:${${prefix}_NAME}>/${output_name}.m
      COMMENT "Copy ${prefix}_NAME documentation file into the output folder"
    )
  endif() # documentation 
  
  # entry point in the mex file + taking care of visibility and symbol clashes.
  if(WIN32)
    set_target_properties(${${prefix}_NAME}
      PROPERTIES 
        DEFINE_SYMBOL "DLL_EXPORT_SYM=__declspec(dllexport)")
  else()
  
    if(HAS_MINUS_PTHREAD AND NOT APPLE)
      # Apparently, compiling with -pthread generated the proper link flags and some defines at compilation
      target_compile_options(${${prefix}_NAME} PRIVATE "-pthread")
    endif()  
  
    if(${prefix}_REDUCE_VISIBILITY)
      # if we do not do that, the symbols linked from eg. boost remain weak and then clash with the ones
      # defined in the matlab process. So by default the symbols are hidden.
      # this also means that for shared libraries (like MEX), the entry point should be explicitely
      # declared with default visibility, otherwise Matlab cannot find the entry point.i
      if(HAS_VISIBILITY_INLINE_HIDDEN)
        target_compile_options(${${prefix}_NAME} PRIVATE "-fvisibility-inlines-hidden")
      endif()
      if(HAS_VISIBILITY_HIDDEN)
        target_compile_options(${${prefix}_NAME} PRIVATE "-fvisibility=hidden")
      endif()    
    endif()
    
    set_target_properties(${${prefix}_NAME} 
      PROPERTIES
        DEFINE_SYMBOL "DLL_EXPORT_SYM=__attribute__ ((visibility (\"default\")))"
    )
  endif()

endfunction(add_matlab_mex)

# this variable will get all Matlab installations found in the current system.
set(_matlab_possible_roots)


# listing the Matlab versions installed on the WIN machine if Matlab_USER_ROOT is not set
if(WIN32)
  
  # for windows
  
  
  if((NOT DEFINED Matlab_USER_ROOT) OR (NOT Matlab_USER_ROOT))
    # if Matlab_USER_ROOT not specified, we look for Matlab installation in the registry
    # if unsuccessful, we look for all known revision and filter the existing ones. 
  
    # testing if we are able to extract the needed information from the registry
    set(matlab_versions_from_registry)
    matlab_extract_all_installed_versions_from_registry(CMAKE_CL_64 matlab_versions_from_registry)
    
    # the returned list is empty, doing the search on all known versions
    if(NOT matlab_versions_from_registry)
      
      if(Matlab_FIND_DEBUG)
        message(STATUS "[Matlab] Search for Matlab from the registry unsuccessful, testing all supported versions")
      endif()
      
      extract_matlab_versions_from_registry_brute_force(matlab_versions_from_registry)
    endif()
        
    # filtering the results with the registry keys
    matlab_get_all_valid_matlab_roots_from_registry(${matlab_versions_from_registry} _matlab_possible_roots)
    
    
    
  elseif(NOT EXISTS ${Matlab_USER_ROOT})
  
    # if Matlab_USER_ROOT specified but erroneous
    if(Matlab_FIND_DEBUG)
      message(WARNING "[Matlab] the specified path for Matlab_USER_ROOT does not exist (${Matlab_USER_ROOT})")
    endif()
  endif()
 
  
else()
  # for linux/osx
  
  if((NOT DEFINED Matlab_USER_ROOT) OR (NOT Matlab_USER_ROOT))

    # if Matlab_USER_ROOT not specified, we look for Matlab from the command line PATH
    # maybe using CMAKE_PROGRAM_PATH to add some more hints
    find_program(
      Matlab_PROGRAM
      "matlab")
  
    if(NOT ${Matlab_PROGRAM})
      #execute_process(COMMAND which matlab OUTPUT_VARIABLE _which_matlab RESULT_VARIABLE _which_matlab_result)
      get_filename_component(Matlab_PROGRAM "matlab" PROGRAM) 
      get_filename_component(Matlab_PROGRAM ${Matlab_PROGRAM} ABSOLUTE) 
      if(EXISTS ${Matlab_PROGRAM})
        if(Matlab_FIND_DEBUG)
          message(STATUS "[Matlab] matlab program result from the command line ${Matlab_PROGRAM}")
        endif()
      else()
        unset(Matlab_PROGRAM)
      endif()

    endif()
  
    if(Matlab_PROGRAM AND EXISTS ${Matlab_PROGRAM})
      if(Matlab_FIND_DEBUG)
        message(STATUS "[Matlab] found from the command line at ${Matlab_PROGRAM}")
      endif()

      # resolve symlinks
      get_filename_component(_matlab_current_location ${Matlab_PROGRAM} REALPATH)
      if(${CMAKE_VERSION} VERSION_LESS "2.8.12")
        set(_directory_alias PATH)
      else()
        set(_directory_alias DIRECTORY)
      endif()
      # get the directory (the command below has to be run twice)
      get_filename_component(_matlab_current_location ${_matlab_current_location} ${_directory_alias})
      get_filename_component(_matlab_current_location ${_matlab_current_location} ${_directory_alias}) # Matlab should be in bin
      list(APPEND _matlab_possible_roots "NOTFOUND" ${_matlab_current_location}) # empty version
    endif()
    
    # on mac, we look for the /Application paths
    # this corresponds to the behaviour on Windows. On Linux, we do not have any other guess.
    if((NOT _matlab_possible_roots) AND APPLE)
      
      matlab_get_supported_releases(_matlab_releases)
      if(Matlab_FIND_DEBUG)
        message(STATUS "[Matlab] Matlab supported versions ${_matlab_releases}. If more version should be supported "
                       "the variable Matlab_ADDITIONAL_VERSIONS can be set according to the documentation")
      endif()

      foreach(_matlab_current_release IN LISTS _matlab_releases)
        set(_matlab_full_string "/Applications/Matlab_${_matlab_current_release}.app")
        if(EXISTS ${_matlab_full_string})
          set(_matlab_current_version)
          matlab_get_version_from_release_name(${_matlab_current_release} _matlab_current_version)
          if(Matlab_FIND_DEBUG)
            message(STATUS "[Matlab] Found version ${_matlab_current_release} (${_matlab_current_version}) in ${_matlab_full_string}")
          endif()
          list(APPEND _matlab_possible_roots ${_matlab_current_version})
          list(APPEND _matlab_possible_roots ${_matlab_full_string})
          unset(_matlab_current_version)
        endif()
        
        unset(_matlab_full_string)
      endforeach(_matlab_current_release)
      unset(_matlab_current_release)
      unset(_matlab_releases) 
    endif()

    # we need to clear Matlab_PROGRAM here
    unset(Matlab_PROGRAM CACHE)
    unset(Matlab_PROGRAM) 
    
  elseif(NOT EXISTS ${Matlab_USER_ROOT})
    # if Matlab_USER_ROOT specified but erroneous
    if(Matlab_FIND_DEBUG)
      message(WARNING "[Matlab] the specified path for Matlab_USER_ROOT does not exist (${Matlab_USER_ROOT})")
    endif()
  else()
    list(APPEND _matlab_possible_roots "NOTFOUND" ${Matlab_USER_ROOT}) # empty version 
  endif()
endif()


if(Matlab_FIND_DEBUG)
  message(STATUS "[Matlab] Matlab root folders are ${_matlab_possible_roots}")
endif()




# take the first possible Matlab root
if(_matlab_possible_roots)
  list(GET _matlab_possible_roots 0 Matlab_VERSION_STRING)
  list(GET _matlab_possible_roots 1 Matlab_ROOT_DIR)
  list(LENGTH _matlab_possible_roots numbers_of_matlab_roots)
  
  # adding a warning in case of ambiguity
  if(numbers_of_matlab_roots GREATER 2 AND Matlab_FIND_DEBUG)
    message(WARNING "[Matlab] Found several distributions of Matlab. Setting the current version to ${Matlab_VERSION_STRING} (located ${Matlab_ROOT_DIR})."
                    " If this is not the desired behaviour, provide the -DMatlab_ROOT on the command line")
  endif()
endif()


if((NOT (DEFINED Matlab_VERSION_STRING)) OR (NOT Matlab_VERSION_STRING))
  set(Matlab_VERSION_STRING "NOTFOUND")
endif()

if(${Matlab_VERSION_STRING} STREQUAL "NOTFOUND")
  if((NOT DEFINED Matlab_PROGRAM) OR (NOT ${Matlab_PROGRAM}) OR (NOT EXISTS ${Matlab_PROGRAM}))
    if(Matlab_FIND_DEBUG)
      message(STATUS "[Matlab] - Unknown version, looking for Matlab under ${Matlab_ROOT_DIR}")
    endif()
    find_program(
      Matlab_PROGRAM
      matlab
      PATHS ${Matlab_ROOT_DIR} ${Matlab_ROOT_DIR}/bin
      DOC "Matlab main program"
      NO_DEFAULT_PATH
    )

     
    
    if(Matlab_PROGRAM)
      set(matlab_list_of_all_versions)
      matlab_get_version_from_matlab_run(${Matlab_PROGRAM} matlab_list_of_all_versions)
    
      list(GET matlab_list_of_all_versions 0 Matlab_VERSION_tmp)
          
      # set the version into the cache
      set(Matlab_VERSION_STRING ${Matlab_VERSION_tmp})# CACHE STRING "Matlab version (automatically determined)")
      list(LENGTH list_of_all_versions list_of_all_versions_length)
      if(${list_of_all_versions_length} GREATER 1)
        message(WARNING "[Matlab] Found several versions, taking the first one (versions found ${list_of_all_versions})")
      endif()
    endif()
  endif()
endif()


if(Matlab_FIND_DEBUG)
  message(STATUS "[Matlab] Current version is ${Matlab_VERSION_STRING} located ${Matlab_ROOT_DIR}")
endif()

if(${Matlab_VERSION_STRING} STREQUAL "NOTFOUND")
  if(Matlab_FIND_DEBUG)
    message(STATUS "Matlab not found. Giving up.")
  endif()
else()

file(TO_CMAKE_PATH ${Matlab_ROOT_DIR} Matlab_ROOT_DIR)

if(CMAKE_SIZEOF_VOID_P EQUAL 4)
  set(_matlab_64Build FALSE)
else()
  set(_matlab_64Build TRUE)
endif()

if(APPLE)
  set(_matlab_bin_prefix "mac") # i should be for intel
  set(_matlab_bin_suffix_32bits "i")
  set(_matlab_bin_suffix_64bits "i64")
elseif(UNIX)
  set(_matlab_bin_prefix "gln")
  set(_matlab_bin_suffix_32bits "x86")
  set(_matlab_bin_suffix_64bits "xa64")
else()
  set(_matlab_bin_prefix "win")
  set(_matlab_bin_suffix_32bits "32")
  set(_matlab_bin_suffix_64bits "64")
endif()



set(Matlab_INCLUDE_DIR_TO_LOOK ${Matlab_ROOT_DIR}/extern/include)
if(_matlab_64Build)
  set(_matlab_current_suffix ${_matlab_bin_suffix_64bits})
else()
  set(_matlab_current_suffix ${_matlab_bin_suffix_32bits})
endif()

set(Matlab_BINARIES_DIR 
    ${Matlab_ROOT_DIR}/bin/${_matlab_bin_prefix}${_matlab_current_suffix} 
    CACHE PATH "Matlab directory for architecture specific binaries" )
set(Matlab_EXTERN_LIBRARY_DIR 
    ${Matlab_ROOT_DIR}/extern/lib/${_matlab_bin_prefix}${_matlab_current_suffix} 
    CACHE PATH "Matlab libraries directory")

if(WIN32)
  set(_matlab_lib_dir_for_search ${Matlab_EXTERN_LIBRARY_DIR}/microsoft)
  set(_matlab_lib_prefix_for_search "lib")
else()
  set(_matlab_lib_dir_for_search ${Matlab_BINARIES_DIR})
  set(_matlab_lib_prefix_for_search "lib")
endif()

unset(_matlab_64Build)


if(NOT DEFINED Matlab_MEX_EXTENSION)
  set(_matlab_mex_extension "")
  matlab_get_mex_suffix(${Matlab_ROOT_DIR} _matlab_mex_extension)

  # This variable goes to the cache.
  set(Matlab_MEX_EXTENSION ${_matlab_mex_extension} CACHE STRING "Extensions for the mex targets (automatically given by Matlab)")
  unset(_matlab_mex_extension)
endif()


if(Matlab_FIND_DEBUG)
  message(STATUS "[Matlab] [DEBUG]_matlab_lib_prefix_for_search = ${_matlab_lib_prefix_for_search} | _matlab_lib_dir_for_search = ${_matlab_lib_dir_for_search}")
endif()

# WARNING: this thing pollutes the CMAKE_FIND_LIBRARY_PREFIXES global variable. 
# Should it be restored afterwards? Is there a more appropriate way to do that?
set(CMAKE_FIND_LIBRARY_PREFIXES ${CMAKE_FIND_LIBRARY_PREFIXES} ${_matlab_lib_prefix_for_search})


set(_matlab_required_variables)


# the MEX library/header are required
    if (Matlab_FIND_DEBUG)
        message(STATUS "[Matlab] Matlab include dir: ${Matlab_INCLUDE_DIRS}")
    endif ()
find_path(
  Matlab_INCLUDE_DIRS
  mex.h
  PATHS ${Matlab_INCLUDE_DIR_TO_LOOK}
  NO_DEFAULT_PATH
  )
list(APPEND _matlab_required_variables Matlab_INCLUDE_DIRS)

    if (Matlab_FIND_DEBUG)
        message(STATUS "[Matlab] Matlab library dir: ${Matlab_MEX_LIBRARY}")
    endif ()
find_library(
  Matlab_MEX_LIBRARY
  mex
  PATHS ${_matlab_lib_dir_for_search}
  NO_DEFAULT_PATH
)
list(APPEND _matlab_required_variables Matlab_MEX_LIBRARY)

# the MEX extension is required
list(APPEND _matlab_required_variables Matlab_MEX_EXTENSION)

# the matlab root is required
list(APPEND _matlab_required_variables Matlab_ROOT_DIR)


# component Mex Compiler
list(FIND Matlab_FIND_COMPONENTS MEX_COMPILER _matlab_find_mex_compiler)
    #if(_matlab_find_mex_compiler GREATER -1)
    if (_matlab_find_mex_compiler EQUAL -1)
  find_program(
    Matlab_MEX_COMPILER
    "mex"
    PATHS ${Matlab_BINARIES_DIR}
    DOC "Matlab MEX compiler"
    NO_DEFAULT_PATH
  )
  
  if(Matlab_MEX_COMPILER)
    set(Matlab_MEX_COMPILER_FOUND TRUE)
  endif()
endif()  
unset(_matlab_find_mex_compiler)

    if (Matlab_FIND_DEBUG)
        message(STATUS "[Matlab] Matlab mex compiler: ${Matlab_MEX_COMPILER}")
    endif ()
# component Matlab program
list(FIND Matlab_FIND_COMPONENTS MAIN_PROGRAM _matlab_find_matlab_program)
    #if(_matlab_find_matlab_program GREATER -1)
    if (_matlab_find_matlab_program EQUAL -1)
  # todo cleanup with code above
  if(NOT DEFINED Matlab_PROGRAM)
    find_program(
      Matlab_PROGRAM
      matlab
      PATHS ${Matlab_ROOT_DIR} ${Matlab_ROOT_DIR}/bin
      DOC "Matlab main program"
      NO_DEFAULT_PATH
    )
  endif()
  if(Matlab_PROGRAM)
    set(Matlab_MAIN_PROGRAM_FOUND TRUE)
  endif()

endif()  
unset(_matlab_find_matlab_program)

    if (Matlab_FIND_DEBUG)
        message(STATUS "[Matlab] Matlab program: ${Matlab_PROGRAM}")
    endif ()
# Component MX library
list(FIND Matlab_FIND_COMPONENTS MX_LIBRARY _matlab_find_mx)
    #if(_matlab_find_mx GREATER -1)
    if (_matlab_find_mx EQUAL -1)
  find_library(
    Matlab_MX_LIBRARY
    mx
    PATHS ${_matlab_lib_dir_for_search}
    NO_DEFAULT_PATH
  )
  
        if (Matlab_FIND_DEBUG)
            message(STATUS "[Matlab] Matlab mx library: ${Matlab_MX_LIBRARY}")
        endif ()
  if(Matlab_MX_LIBRARY)
    set(Matlab_MX_LIBRARY_FOUND TRUE)
  endif()
endif()
unset(_matlab_find_mx)

    if (Matlab_FIND_DEBUG)
        message(STATUS "[Matlab] Matlab mx library: ${_matlab_lib_dir_for_search}")
    endif ()

# Component ENG library
list(FIND Matlab_FIND_COMPONENTS ENG_LIBRARY _matlab_find_eng)
    #if(_matlab_find_eng GREATER -1)
    if (_matlab_find_eng EQUAL -1)
  find_library(
    Matlab_ENG_LIBRARY
    eng
    PATHS ${_matlab_lib_dir_for_search}
    NO_DEFAULT_PATH
  )
  if(Matlab_ENG_LIBRARY)
    set(Matlab_ENG_LIBRARY_FOUND TRUE)
  endif()
endif()
unset(_matlab_find_eng)


unset(_matlab_lib_dir_for_search)
endif(${Matlab_VERSION_STRING} STREQUAL "NOTFOUND")

set(Matlab_LIBRARIES ${Matlab_MEX_LIBRARY} ${Matlab_MX_LIBRARY} ${Matlab_ENG_LIBRARY})

if (Matlab_FIND_DEBUG)
    message(STATUS "[Matlab] set matlab libraries: ${Matlab_LIBRARIES}")
endif ()
unset(_matlab_required_variables)
unset(_matlab_bin_prefix)
unset(_matlab_bin_suffix_32bits)
unset(_matlab_bin_suffix_64bits)
unset(_matlab_current_suffix)
unset(_matlab_lib_dir_for_search)
unset(_matlab_lib_prefix_for_search)

if (Matlab_FIND_DEBUG)
    message(STATUS "[Matlab] matlab found: ${Matlab_INCLUDE_DIRS} and ${Matlab_LIBRARIES}")
endif ()
if(Matlab_INCLUDE_DIRS AND Matlab_LIBRARIES)
    # handle REQUIRED and QUIET options
    include(FindPackageHandleStandardArgs)
    if (CMAKE_VERSION LESS 2.8.3)
        find_package_handle_standard_args(MATLAB DEFAULT_MSG Matlab_INCLUDE_DIRS Matlab_LIBRARIES Matlab_ROOT_DIR)
    else ()
        find_package_handle_standard_args(MATLAB REQUIRED_VARS Matlab_INCLUDE_DIRS Matlab_LIBRARIES VERSION_VAR Matlab_ROOT_DIR)
    endif ()


  mark_as_advanced(
    Matlab_LIBRARIES
    Matlab_MEX_LIBRARY
    Matlab_MX_LIBRARY
    Matlab_ENG_LIBRARY
    Matlab_INCLUDE_DIRS
    Matlab_FOUND
    Matlab_USER_ROOT
    Matlab_ROOT_DIR
    Matlab_VERSION_STRING
    Matlab_PROGRAM
    Matlab_MEX_EXTENSION
    Matlab_BINARIES_DIR
  )
endif()




