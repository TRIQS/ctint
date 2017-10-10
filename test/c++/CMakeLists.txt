# Copy h5 files to binary dir
#file(GLOB all_h5_files *.h5)
#file(COPY ${all_h5_files} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} FILES_MATCHING PATTERN "*.h5")

# NFFT Test
add_executable(nfft ${CMAKE_CURRENT_SOURCE_DIR}/nfft.cpp)
target_link_libraries (nfft ctint_c gtest ${NFFT_LIBRARIES})
add_test(nfft ${CMAKE_CURRENT_BINARY_DIR}/nfft)

## Anderson Test
#add_executable(anderson ${CMAKE_CURRENT_SOURCE_DIR}/anderson.cpp)
#target_link_libraries (anderson ctint_c ${TRIQS_LIBRARY_ALL} ${TRIQS_LIBRARY_GTEST})
#add_test(anderson ${CMAKE_CURRENT_BINARY_DIR}/anderson)
