# Version number
set (CTINT_VERSION "1.5")

# Start configuration
cmake_minimum_required(VERSION 2.8.12)
project(triqs_ctint CXX)

# Default to Release build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Type of build" FORCE)
endif()
message( STATUS "-------- BUILD-TYPE: ${CMAKE_BUILD_TYPE} --------")

# Options for the compilation ...
option(GTAU_IS_COMPLEX "If ON, the non-interacting Green-function G_0(tau) is complex" OFF)
option(INTERACTION_IS_COMPLEX "If ON, the interaction vertices (dynamic and static) are complex" OFF)

# Build static libraries
set(BUILD_SHARED_LIBS OFF)

# Global options and definitions
add_definitions(-Wall)
add_compile_options($<$<CONFIG:Release>:-march=native>
		    $<$<CONFIG:Debug>:-Og>)

# Load Dependencies
find_package(TRIQS ${CTINT_VERSION} EXACT REQUIRED)
find_package(Cpp2Py ${CTINT_VERSION} EXACT REQUIRED)
find_package(NFFT REQUIRED)

# Default Install directory to TRIQS_ROOT if not given or invalid.
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT OR (NOT IS_ABSOLUTE ${CMAKE_INSTALL_PREFIX}))
  message(STATUS "No install prefix given (or invalid). Defaulting to TRIQS_ROOT")
  set(CMAKE_INSTALL_PREFIX ${TRIQS_ROOT} CACHE PATH "default install path" FORCE)
endif()
message(STATUS "-------- INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX} --------")

# Fetch the git hash of ctint
triqs_get_git_hash_of_source_dir(CTINT_GIT_HASH)
message(STATUS "Git hash: ${CTINT_GIT_HASH}")

# Build and install the app4triqs library
add_subdirectory(c++/triqs_ctint)

# Build and install the ctint python module
if(${TRIQS_WITH_PYTHON_SUPPORT})
 add_subdirectory(python/triqs_ctint)
endif()

# Tests
enable_testing()
add_subdirectory(test)

# Additional configuration files
add_subdirectory(share)

option(Build_Documentation "Build documentation" OFF)
if(${Build_Documentation})
  if(NOT ${TRIQS_WITH_DOCUMENTATION})
    message("Error: TRIQS library has not been compiled with its documentation")
  else()
    add_subdirectory(doc)
  endif()
endif()
