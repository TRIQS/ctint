# Copy h5 files to binary dir
file(GLOB all_h5_files *.h5)
file(COPY ${all_h5_files} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

set(all_tests util gf_real_imag nfft_buf)

foreach(t ${all_tests})
  add_executable(${t} ${CMAKE_CURRENT_SOURCE_DIR}/${t}.cpp)
  target_link_libraries(${t} ctint_c gtest nfft)
  add_test(${t} ${CMAKE_CURRENT_BINARY_DIR}/${t})
endforeach()
