# Version number of CTINT
set (CTINT_VERSION "1.5")

# Start configuration
cmake_minimum_required(VERSION 2.8)
project(triqs_ctint CXX)

# Default to Release build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Type of build" FORCE)
endif()
message( STATUS "-------- BUILD-TYPE: ${CMAKE_BUILD_TYPE} -------------")

option(GTAU_IS_COMPLEX "If ON, the non-interacting Green-function G_0(tau) is complex" OFF)
option(INTERACTION_IS_COMPLEX "If ON, the interaction vertices (dynamic and static) are complex" OFF)

# Use shared libraries
set(BUILD_SHARED_LIBS ON)

# Load TRIQS and Cpp2Py
find_package(TRIQS REQUIRED)
find_package(Cpp2Py REQUIRED)

# We will need the FindNNFT of triqs, so we adjust CMake Module Path to find it
list(APPEND CMAKE_MODULE_PATH ${TRIQS_ROOT}/share/cmake)

# Check that ctint is compatible with the TRIQS version
if(NOT CTINT_VERSION EQUAL TRIQS_VERSION)
  message(FATAL_ERROR "The application version is not compatible with the TRIQS library (TRIQS library version: ${TRIQS_VERSION} while this application version: ${CTINT_VERSION})")
endif()

# Default Install directory to TRIQS_ROOT
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT OR (NOT IS_ABSOLUTE ${CMAKE_INSTALL_PREFIX}))
  message(STATUS " No install prefix given (or invalid). Defaulting to TRIQS_ROOT")
  set(CMAKE_INSTALL_PREFIX ${TRIQS_ROOT} CACHE PATH "default install path" FORCE)
endif()
message(STATUS "-------- CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX} -------------")

triqs_get_git_hash_of_source_dir(CTINT_GIT_HASH)

add_compile_options(    -std=c++17 -DHAS_OPTIONAL_HEADER 
  			-DCTINT_GIT_HASH=${CTINT_GIT_HASH} -DTRIQS_GIT_HASH=${TRIQS_GIT_HASH}
			$<$<CONFIG:Release>:-march=native>
			$<$<CONFIG:Debug>:-DDEBUG_CTINT>
			$<$<BOOL:GTAU_IS_COMPLEX>:-DGTAU_IS_COMPLEX>
			$<$<BOOL:INTERACTION_IS_COMPLEX>:-DINTERACTION_IS_COMPLEX>
			)

# Common to all targets (ctint, python, tests).
include_directories(c++/triqs_ctint)

# To be fixed!!
triqs_set_rpath_for_all_targets()

add_subdirectory(c++/triqs_ctint)
add_subdirectory(python/triqs_ctint)

# Tests
enable_testing()
add_subdirectory(test)

# Build the documentation
option(Build_Documentation "Build documentation" OFF)
if(${Build_Documentation})
  if(NOT TRIQS_WITH_DOCUMENTATION)
    message("Error: TRIQS library has not been compiled with its documentation")
  endif()
  add_subdirectory(doc)
endif()
