FILE(GLOB_RECURSE cpp_sources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)
add_library(ctint_c ${cpp_sources})
target_link_libraries(ctint_c PUBLIC triqs)

target_link_libraries(ctint_c PRIVATE nfft)

# Install
install(TARGETS ctint_c DESTINATION lib)
