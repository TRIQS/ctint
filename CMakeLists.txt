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
add_compile_options(-Wall
		     $<$<CONFIG:Debug>:-Og>
		     $<$<CONFIG:Debug>:-ggdb3>
		     )
# Load Dependencies
find_package(TRIQS 2.0 REQUIRED)
find_package(Cpp2Py 1.5 REQUIRED)
find_package(NFFT REQUIRED)
install(TARGETS nfft EXPORT ctint-targets)

# Default Install directory to TRIQS_ROOT if not given or invalid.
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT OR (NOT IS_ABSOLUTE ${CMAKE_INSTALL_PREFIX}))
  message(STATUS "No install prefix given (or invalid). Defaulting to TRIQS_ROOT")
  set(CMAKE_INSTALL_PREFIX ${TRIQS_ROOT} CACHE PATH "default install path" FORCE)
endif()
message(STATUS "-------- INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX} --------")

# Use git to detect CTINT version and git hash
execute_process(COMMAND bash "-c" "${TRIQS_GIT_EXECUTABLE} describe | sed \"s/-\\([0-9]*\\)-.*/.\\1/g\""
 WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
 OUTPUT_VARIABLE   CTINT_VERSION
 OUTPUT_STRIP_TRAILING_WHITESPACE
)
if(NOT CTINT_VERSION)
  set(CTINT_VERSION 0.8)
endif()
triqs_get_git_hash_of_source_dir(CTINT_GIT_HASH)
MESSAGE(STATUS "ctint version : ${CTINT_VERSION}")
message(STATUS "Git hash: ${CTINT_GIT_HASH}")

# Build and install the ctint library
add_subdirectory(c++/triqs_ctint)

# Build and install the ctint python module
if(${TRIQS_WITH_PYTHON_SUPPORT})
 add_subdirectory(python/triqs_ctint)
endif()

# Tests
option(Build_Tests "Build tests" ON)
if(${Build_Tests})
  enable_testing()
  add_subdirectory(test)
endif()

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

option(BUILD_DEBIAN_PACKAGE "Build a deb package" OFF)
if(BUILD_DEBIAN_PACKAGE)
  if(NOT CMAKE_INSTALL_PREFIX STREQUAL "/usr")
    message(FATAL_ERROR "CMAKE_INSTALL_PREFIX must be /usr for packaging")
  endif()
  SET(CPACK_GENERATOR "DEB")
  SET(CPACK_PACKAGE_VERSION ${CTINT_VERSION})
  SET(CPACK_PACKAGE_CONTACT "https://github.com/TRIQS/ctint")
  EXECUTE_PROCESS(COMMAND dpkg --print-architecture OUTPUT_VARIABLE CMAKE_DEBIAN_PACKAGE_ARCHITECTURE OUTPUT_STRIP_TRAILING_WHITESPACE)
  SET(CPACK_DEBIAN_PACKAGE_DEPENDS "triqs (>= 2.0), cpp2py (>= 1.5)")
  SET(CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)
  SET(CPACK_DEBIAN_PACKAGE_GENERATE_SHLIBS ON)
  INCLUDE(CPack)
endif()
