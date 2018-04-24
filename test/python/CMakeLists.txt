# If ASAN=ON use minimal runtime of UBSAN if it exists
if(${ASAN} AND UBSAN_MINIMAL_RT_LIBRARY)
  set(UBSAN_RT_LIBRARY ${UBSAN_MINIMAL_RT_LIBRARY})
endif()

add_subdirectory(static)
add_subdirectory(dynamic)
