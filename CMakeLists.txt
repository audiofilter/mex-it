cmake_minimum_required(VERSION 2.8.4)

if (POLICY CMP0054)
  cmake_policy(SET CMP0054 NEW)
endif()

if (POLICY CMP0042)
  cmake_policy(SET CMP0042 NEW)
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake_modules")

add_subdirectory(test)