# Copy h5 files to binary dir
file(GLOB all_h5_files *.h5)
file(COPY ${all_h5_files} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# NFFT Test
add_executable(test_nfft ${CMAKE_CURRENT_SOURCE_DIR}/nfft.cpp)
target_link_libraries (test_nfft triqs_ctint_c gtest ${NFFT_LIBRARIES})
add_test(test_nfft ${CMAKE_CURRENT_BINARY_DIR}/test_nfft)

set(all_tests util gf_real_imag)

foreach(t ${all_tests})
  add_executable(${t} ${CMAKE_CURRENT_SOURCE_DIR}/${t}.cpp)
  target_link_libraries(${t} triqs_ctint_c gtest)
  add_test(${t} ${CMAKE_CURRENT_BINARY_DIR}/${t})
endforeach()
