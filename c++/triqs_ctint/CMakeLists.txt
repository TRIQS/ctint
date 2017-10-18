FILE(GLOB_RECURSE cpp_sources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)
add_library(triqs_ctint_c ${cpp_sources})
target_link_libraries(triqs_ctint_c PUBLIC triqs)

# Find the NFFT Library
message( STATUS "-------- NFFT detection -------------")
find_package(NFFT REQUIRED)
target_include_directories(triqs_ctint_c SYSTEM PRIVATE ${NFFT_INCLUDE_DIR})
if(NFFT_OLD_API)
  target_compile_definitions(triqs_ctint_c PUBLIC NFFT_OLD_API)
endif(NFFT_OLD_API)
target_link_libraries(triqs_ctint_c PUBLIC ${NFFT_LIBRARIES})

# Install
install(TARGETS triqs_ctint_c DESTINATION lib)
