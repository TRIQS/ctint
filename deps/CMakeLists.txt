include(external_dependency.cmake)

# Add your dependencies with the function
#
#   external_dependency(name
#       [VERSION <version-number>]
#       [GIT_REPO <url>]
#       [GIT_TAG <tag>]
#       [BUILD_ALWAYS]
#       [EXCLUDE_FROM_ALL]
#   )
#
# Resolve the dependency using the following steps in order.
# If a step was successful, skip the remaining ones.
#
#  1. Use find_package(name [<version-number>])
#     to locate the package in the system.
#     Skip this step if Build_Deps option is set.
#  2. Try to find a directory containing the sources
#     at ${CMAKE_CURRENT_SOURCE_DIR}/name and
#     ${CMAKE_SOURCE_DIR}/deps/name. If found
#     build it as a cmake sub-project.
#  3. If GIT_REPO is provided, git clone the sources,
#     and build them as a cmake sub-project.
#
# Addtional options:
#
#  GIT_TAG - Use this keyword to specify the git-tag, branch or commit hash
#
#  BUILD_ALWAYS - If set, this dependency will always be built from source
#                 and will never be searched in the system.
#
#  EXCLUDE_FROM_ALL - If set, targets of the dependency cmake subproject
#                     will not be included in the ALL target of the project.
#                     In particular the dependency will not be installed.

if(NOT DEFINED Build_Deps)
  set(Build_Deps "Always" CACHE STRING "Do we build dependencies from source? [Never/Always/IfNotFound]")
else()
  set(Build_Deps_Opts "Never" "Always" "IfNotFound")
  if(NOT ${Build_Deps} IN_LIST Build_Deps_Opts)
    message(FATAL_ERROR "Build_Deps option should be either 'Never', 'Always' or 'IfNotFound'")
  endif()
  set(Build_Deps ${Build_Deps} CACHE STRING "Do we build dependencies from source? [Never/Always/IfNotFound]")
  if(NOT IS_SUBPROJECT AND NOT Build_Deps STREQUAL "Always" AND (ASAN OR UBSAN))
    message(WARNING "For builds with llvm sanitizers (ASAN/UBSAN) it is recommended to use -DBuild_Deps=Always to avoid false positives.")
  endif()
endif()

# -- triqs_hartree_fock --
if(PythonSupport)
  external_dependency(triqs_hartree_fock
    GIT_REPO https://github.com/TRIQS/hartree_fock
    VERSION 4.0
    GIT_TAG unstable
  )
endif()

# -- finufft --
# BUILD_ALWAYS: c++/triqs_ctint absorbs finufft's objects into libtriqs_ctint_c.a via
# $<TARGET_OBJECTS:>, which does not work on the IMPORTED target a system finufft would give us.
# Note: finufft fetches xsimd via CPM, which we use for ctint as well.
set(FINUFFT_USE_OPENMP OFF)
# Default to a portable finufft build: its -march=native default crashes (SIGILL) when a cross-CPU build cache is reused. Override with -DFINUFFT_ARCH_FLAGS=native for a local perf build.
set(FINUFFT_ARCH_FLAGS "" CACHE STRING "Architecture flags for finufft (empty = portable baseline)")
set(FINUFFT_USE_DUCC0 ON)
set(FINUFFT_INTERPROCEDURAL_OPTIMIZATION ON)
external_dependency(finufft
  GIT_REPO https://github.com/diamondinoia/finufft
  GIT_TAG experimental
  BUILD_ALWAYS
  EXCLUDE_FROM_ALL
)

# -- polyfit (also provides poet as transitive dependency) --
set(POLYFIT_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(POLYFIT_BUILD_BENCHMARKS OFF CACHE BOOL "" FORCE)
external_dependency(polyfit
  GIT_REPO https://github.com/DiamonDinoia/polyfit
  GIT_TAG main
  BUILD_ALWAYS
  EXCLUDE_FROM_ALL
)


# -- GTest --
external_dependency(GTest
  GIT_REPO https://github.com/google/googletest
  GIT_TAG main
  BUILD_ALWAYS
  EXCLUDE_FROM_ALL
)

# -- google-bench --
if(NOT IS_SUBPROJECT AND Build_Benchs AND NOT ASAN AND NOT UBSAN)
  set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "")
  external_dependency(google-bench
    GIT_REPO https://github.com/google/benchmark
    GIT_TAG main
    BUILD_ALWAYS
    EXCLUDE_FROM_ALL
  )
endif()
