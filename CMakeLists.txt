# Version number of the application
set (CTINT_VERSION "1.5")
set (CTINT_RELEASE "1.5.0")

# Append triqs installed files to the cmake load path
list(APPEND CMAKE_MODULE_PATH ${TRIQS_PATH}/share/triqs/cmake)

# start configuration 
cmake_minimum_required(VERSION 2.8)
project(ctint CXX)
set(CMAKE_BUILD_TYPE Release)
enable_testing()

# We use shared libraries
option(BUILD_SHARED_LIBS "Build shared libraries" ON)

# Find NFFT Library
message( STATUS "-------- NFFT detection -------------")
find_package(NFFT REQUIRED)
include_directories(SYSTEM ${NFFT_INCLUDE_DIR})
if(NFFT_OLD_API)
  add_definitions(-DNFFT_OLD_API)
endif(NFFT_OLD_API)

# Load TRIQS, including all predefined variables from TRIQS installation
message( STATUS "-------- TRIQS detection ------------")
find_package(TRIQS REQUIRED)
include_directories(${TRIQS_INCLUDE_ALL})
add_definitions(${TRIQS_CXX_DEFINITIONS}) 

# Check that versions are compatible
if(NOT CTINT_VERSION EQUAL TRIQS_VERSION)
  message(FATAL_ERROR "The application version is not compatible with the TRIQS library (TRIQS library version: ${TRIQS_VERSION} while this application version: ${CTHYB_SEGMENT_VERSION})")
endif()

option(DEBUG_CTINT "Debug mode" OFF)
option(VERBOSE_CTINT "Verbose mode" OFF)

# We want to be installed in the TRIQS tree
set(CMAKE_INSTALL_PREFIX ${TRIQS_PATH})

# Compile C++ code
message( STATUS "-------- Compiling C++ Code ---------")
include_directories(c++)
add_subdirectory(c++)

# Python interface
if (${TRIQS_WITH_PYTHON_SUPPORT})
  message( STATUS "---------- Python Interface ---------")
  add_subdirectory(python)
  option(DocWithCpp2doc "Use c++2doc to generate C++ documentation (requires libclang)" OFF)
  option(BUILD_DOC "Build documentation" OFF)
  if(${BUILD_DOC})
    if(NOT TRIQS_WITH_DOCUMENTATION)
      message("Error: TRIQS library has not been compiled with its documentation")
    endif()
    add_subdirectory(doc)
  endif(${BUILD_DOC})
endif()

# Add tests
add_subdirectory(test)
