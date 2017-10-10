# Version number of CTINT
set (CTINT_VERSION "1.5")

# start configuration
cmake_minimum_required(VERSION 2.8)
project(ctint CXX)

# Default to Release build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Type of build" FORCE)
endif()
message( STATUS "-------- BUILD-TYPE: ${CMAKE_BUILD_TYPE} -------------")

# Use shared libraries
set(BUILD_SHARED_LIBS ON)

# Load TRIQS, (and Cpp2Py automatically)
find_package(TRIQS REQUIRED)

# We will need the FindNNFT of triqs, so we adjust CMake Module Path to find it
list(APPEND CMAKE_MODULE_PATH ${TRIQS_ROOT}/share/cmake)

# Check that ctint is compatible with the TRIQS version
if(NOT CTINT_VERSION EQUAL TRIQS_VERSION)
  message(FATAL_ERROR "The application version is not compatible with the TRIQS library (TRIQS library version: ${TRIQS_VERSION} while this application version: ${CTINT_VERSION})")
endif()

get_git_hash_of_source_dir(GIT_HASH)

add_compile_options(    -std=c++17 -DHAS_OPTIONAL_HEADER -DGIT_SHA1=${GIT_HASH}
			"$<$<CONFIG:Release>:-march=native>"
			"$<$<CONFIG:Debug>:-DDEBUG_CTINT>"
			)

# Common to all targets (ctint, python, tests).
include_directories(c++)

# To be fixed!!
triqs_set_rpath_for_all_targets()

add_subdirectory(c++)
add_subdirectory(python)

# Tests
enable_testing()
add_subdirectory(test)

# Build the documentation
option(BUILD_DOC "Build documentation" OFF)
if(${BUILD_DOC})
  if(NOT TRIQS_WITH_DOCUMENTATION)
    message("Error: TRIQS library has not been compiled with its documentation")
  endif()
  add_subdirectory(doc)
endif()
