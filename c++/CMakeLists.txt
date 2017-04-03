if (${DEBUG_CTINT})
  message(STATUS "DEBUG MODE")
  add_definitions(-DDEBUG_CTINT)
endif()
if (${VERBOSE_CTINT})
  message(STATUS "VERBOSE MODE")
  add_definitions(-DVERBOSE_CTINT)
endif()

FILE(GLOB_RECURSE cpp_sources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)
#add_library(ctint_c ${cpp_sources})
add_library(ctint_c ${cpp_sources})
target_link_libraries(ctint_c ${TRIQS_LIBRARY_ALL} ${NFFT_LIBRARIES})
triqs_set_rpath_for_target(ctint_c)

install(TARGETS ctint_c DESTINATION lib)

