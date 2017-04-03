
# Copy h5 files to binary dir
FILE(GLOB all_h5_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.h5)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/${all_h5_files} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
triqs_add_python_test(densdens)
triqs_add_python_test(jperp)
triqs_add_python_test(all)
triqs_add_python_test(all_continuousBoson)

