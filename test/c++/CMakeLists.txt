# Copy h5 files to binary dir
FILE(GLOB all_h5_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.h5)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/${all_h5_files} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

include_directories(${CMAKE_SOURCE_DIR}/c++)

# NFFT Test
add_executable(nfft ${CMAKE_CURRENT_SOURCE_DIR}/nfft.cpp)
target_link_libraries (nfft ${TRIQS_LIBRARY_ALL} ${TRIQS_LIBRARY_GTEST} ${NFFT_LIBRARIES})
add_test(nfft ${CMAKE_CURRENT_BINARY_DIR}/nfft)

## Anderson Test
#add_executable(anderson ${CMAKE_CURRENT_SOURCE_DIR}/anderson.cpp)
#target_link_libraries (anderson ctint_c ${TRIQS_LIBRARY_ALL} ${TRIQS_LIBRARY_GTEST})
#add_test(anderson ${CMAKE_CURRENT_BINARY_DIR}/anderson)

#triqs_add_all_cpp_tests_in_current_dir()
