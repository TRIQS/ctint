if (CMAKE_BUILD_TYPE MATCHES Debug)
  add_definitions(-DDEBUG_CTINT)
endif()

FILE(GLOB_RECURSE cpp_sources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)
add_library(ctint_c ${cpp_sources})
target_link_libraries(ctint_c ${TRIQS_LIBRARY_ALL} ${NFFT_LIBRARIES})
triqs_set_rpath_for_target(ctint_c)

install(TARGETS ctint_c DESTINATION lib)

