# Version number of CTINT
set (CTINT_VERSION "1.5")

# start configuration
cmake_minimum_required(VERSION 2.8)
project(ctint CXX)

# Adjust CMake Module Path and includes
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake ${TRIQS_ROOT}/share/cmake)
include(GetGitRevisionDescription)

# Default to Release build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Type of build" FORCE)
endif()
message( STATUS "-------- BUILD-TYPE: ${CMAKE_BUILD_TYPE} -------------")

enable_testing()

# Use shared libraries
set(BUILD_SHARED_LIBS ON)

# Load TRIQS, (and Cpp2Py automatically)
find_package(TRIQS REQUIRED)

# Check that ctint is compatible with the TRIQS version
if(NOT CTINT_VERSION EQUAL TRIQS_VERSION)
  message(FATAL_ERROR "The application version is not compatible with the TRIQS library (TRIQS library version: ${TRIQS_VERSION} while this application version: ${CTINT_VERSION})")
endif()

add_compile_options(    -std=c++17 -DHAS_OPTIONAL_HEADER -DGIT_SHA1=${GIT_SHA1}
			"$<$<CONFIG:Release>:-march=native>"
			"$<$<CONFIG:Debug>:-DDEBUG_CTINT>"
			)

# Sets the rpath for all targets, as for all TRIQS applications.
triqs_set_rpath_for_all_targets()

add_subdirectory(c++)
add_subdirectory(python)
add_subdirectory(test)

include_directories(c++)

# Build the documentation
option(BUILD_DOC "Build documentation" OFF)
if(${BUILD_DOC})
  if(NOT TRIQS_WITH_DOCUMENTATION)
    message("Error: TRIQS library has not been compiled with its documentation")
  endif()
  add_subdirectory(doc)
endif()
