file(GLOB_RECURSE sources *.cpp)
add_library(ctint_c ${sources})

# Link against triqs and use headers with SYSTEM flag for better warning messages
target_link_libraries(ctint_c PUBLIC triqs)
target_include_directories(ctint_c SYSTEM PUBLIC ${TRIQS_ROOT}/include)

target_link_libraries(ctint_c PRIVATE nfft)

target_compile_options(ctint_c PUBLIC -std=c++17 -fPIC)
target_include_directories(ctint_c PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/c++>)
target_compile_definitions(ctint_c PUBLIC
			    CTINT_GIT_HASH=${CTINT_GIT_HASH}
			    TRIQS_GIT_HASH=${TRIQS_GIT_HASH}
			    $<$<CONFIG:Debug>:DEBUG_CTINT>
			    $<$<CONFIG:Debug>:TRIQS_ARRAYS_ENFORCE_BOUNDCHECK>
			    $<$<BOOL:${GTAU_IS_COMPLEX}>:GTAU_IS_COMPLEX>
			    $<$<BOOL:${INTERACTION_IS_COMPLEX}>:INTERACTION_IS_COMPLEX>
			  )

# Install library and headers
install(TARGETS ctint_c EXPORT ctint-targets DESTINATION lib)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} DESTINATION include FILES_MATCHING PATTERN "*.hpp" PATTERN "*.h")

# ========= Static Analyzer Checks ==========

# Locate static analyzer tools
find_program(CPPCHECK_EXECUTABLE NAMES "cppcheck" PATHS ENV PATH)
find_program(CLANG_TIDY_EXECUTABLE NAMES "clang-tidy" PATHS ENV PATH)

# Run clang-tidy if found
if(CLANG_TIDY_EXECUTABLE)
  message(STATUS "clang-tidy found: ${CLANG_TIDY_EXECUTABLE}")
  set_target_properties(ctint_c PROPERTIES CXX_CLANG_TIDY "${CLANG_TIDY_EXECUTABLE}")
else()
  message(STATUS "clang-tidy not found in $PATH. Please consider installing clang-tidy for additional checks!")
endif()

# Run cppcheck if found
if(CPPCHECK_EXECUTABLE)
  message(STATUS "cppcheck found: ${CPPCHECK_EXECUTABLE}")
  add_custom_command(
	TARGET ctint_c
	COMMAND ${CPPCHECK_EXECUTABLE}
	--enable=warning,style,performance,portability
        --std=c++14
	--template=gcc
        --verbose
        --quiet
        ${sources}
  )
else()
  message(STATUS "cppcheck not found in $PATH. Please consider installing cppcheck for additional checks!")
endif()


# ========= Dynamic Analyzer Checks ==========

# Address Sanitizer
option(ASAN OFF "Compile library and executables with LLVM Address Sanitizer")
if(ASAN)
  list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/share/cmake)
  find_package(libasan_rt REQUIRED)
  target_compile_options(ctint_c PUBLIC -fsanitize=address -fno-omit-frame-pointer)
  target_link_libraries(ctint_c INTERFACE "-fsanitize=address -fno-omit-frame-pointer")
endif()

# Undefined Behavior Sanitizer
option(UBSAN OFF "Compile library and executables with LLVM Undefined Behavior Sanitizer")
if(UBSAN)
  list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/share/cmake)
  find_package(libubsan_rt REQUIRED)
  target_compile_options(ctint_c PUBLIC -fsanitize=undefined -fno-omit-frame-pointer -fno-sanitize=vptr)
  target_link_libraries(ctint_c INTERFACE "-fsanitize=undefined -fno-omit-frame-pointer")
endif()
