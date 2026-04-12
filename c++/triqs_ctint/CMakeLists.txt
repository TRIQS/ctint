file(GLOB_RECURSE sources *.cpp)

# ========= Variant Library Function ==========

function(add_ctint_library TARGET_NAME)
  cmake_parse_arguments(ARG "" "" "COMPILE_DEFS" ${ARGN})

  add_library(${TARGET_NAME} ${sources})
  add_library(${PROJECT_NAME}::${TARGET_NAME} ALIAS ${TARGET_NAME})

  # Link against triqs and enable warnings
  # finufft, xsimd, poet come transitively through triqs (TRIQS_WITH_NFFT)
  target_link_libraries(${TARGET_NAME}
    PUBLIC
      triqs
      $<BUILD_INTERFACE:${PROJECT_NAME}_performance_flags>
    PRIVATE
      $<BUILD_INTERFACE:${PROJECT_NAME}_warnings>
  )

  # Configure target and compilation
  set_target_properties(${TARGET_NAME} PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    VERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
  )
  target_include_directories(${TARGET_NAME} PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/c++>)
  target_include_directories(${TARGET_NAME} SYSTEM INTERFACE $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include>)
  target_link_directories(${TARGET_NAME} PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/c++>)
  target_link_directories(${TARGET_NAME} INTERFACE $<INSTALL_INTERFACE:${CMAKE_INSTALL_FULL_LIBDIR}>)
  target_compile_definitions(${TARGET_NAME} PUBLIC
    CTINT_GIT_HASH=${PROJECT_GIT_HASH}
    TRIQS_GIT_HASH=${TRIQS_GIT_HASH}
    $<$<CONFIG:Debug>:CTINT_DEBUG>
    $<$<CONFIG:Debug>:TRIQS_DEBUG>
    $<$<CONFIG:Debug>:TRIQS_ARRAYS_ENFORCE_BOUNDCHECK>
    ${ARG_COMPILE_DEFS}
  )

  # Sanitizers
  foreach(san IN ITEMS asan ubsan msan tsan)
    string(TOUPPER ${san} SAN)
    if(${SAN})
      if(NOT TARGET ${san})
        find_package(sanitizer REQUIRED COMPONENTS ${san})
      endif()
      target_link_libraries(${TARGET_NAME} PUBLIC $<BUILD_INTERFACE:${san}>)
    endif()
  endforeach()

  install(TARGETS ${TARGET_NAME} EXPORT ${PROJECT_NAME}-targets DESTINATION ${CMAKE_INSTALL_LIBDIR})
endfunction()

# ========= Dynamic Analyzer Options ==========

option(ASAN OFF "Compile library and executables with LLVM Address Sanitizer")
option(UBSAN OFF "Compile library and executables with LLVM Undefined Behavior Sanitizer")
option(MSAN OFF "Compile library and executables with LLVM Memory Sanitizer")
option(TSAN OFF "Compile library and executables with LLVM Thread Sanitizer")

# ========= Build Library Variants ==========

# Always build the real variant
add_ctint_library(${PROJECT_NAME}_c)

# Optionally build complex G(tau) variant
if(GTAU_IS_COMPLEX)
  add_ctint_library(${PROJECT_NAME}_c_complex_gtau COMPILE_DEFS GTAU_IS_COMPLEX)
endif()

# Optionally build complex G(tau) + interaction variant
if(INTERACTION_IS_COMPLEX)
  if(NOT GTAU_IS_COMPLEX)
    message(FATAL_ERROR "INTERACTION_IS_COMPLEX requires GTAU_IS_COMPLEX=ON")
  endif()
  add_ctint_library(${PROJECT_NAME}_c_complex_all COMPILE_DEFS GTAU_IS_COMPLEX INTERACTION_IS_COMPLEX)
endif()

# Install headers (once)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} DESTINATION include FILES_MATCHING PATTERN "*.hpp" PATTERN "*.h")

# finufft, xsimd, poet are now provided by triqs — no need to install separately

# ========= Static Analyzer Checks ==========

option(ANALYZE_SOURCES OFF "Run static analyzer checks if found (clang-tidy, cppcheck)")
if(ANALYZE_SOURCES)

  # Locate static analyzer tools
  find_program(CPPCHECK_EXECUTABLE NAMES "cppcheck" PATHS ENV PATH)
  find_program(CLANG_TIDY_EXECUTABLE NAMES "clang-tidy" PATHS ENV PATH)

  # Run clang-tidy if found
  if(CLANG_TIDY_EXECUTABLE)
    message(STATUS "clang-tidy found: ${CLANG_TIDY_EXECUTABLE}")
    set_target_properties(${PROJECT_NAME}_c PROPERTIES CXX_CLANG_TIDY "${CLANG_TIDY_EXECUTABLE}")
  else()
    message(STATUS "clang-tidy not found in $PATH. Please consider installing clang-tidy for additional checks!")
  endif()

  # Run cppcheck if found
  if(CPPCHECK_EXECUTABLE)
    message(STATUS "cppcheck found: ${CPPCHECK_EXECUTABLE}")
    add_custom_command(
	TARGET ${PROJECT_NAME}_c
	COMMAND ${CPPCHECK_EXECUTABLE}
	  --enable=warning,style,performance,portability
	  --std=c++23
	  --template=gcc
	  --verbose
	  --force
	  --quiet
          ${sources}
	WORKING_DIRECTORY
	  ${CMAKE_CURRENT_SOURCE_DIR}
    )
  else()
    message(STATUS "cppcheck not found in $PATH. Please consider installing cppcheck for additional checks!")
  endif()

endif()
