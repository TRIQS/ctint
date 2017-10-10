# Copy h5 files to binary dir
#FILE(GLOB all_h5_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.h5)
#file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/${all_h5_files} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} FILES_MATCHING PATTERN "*.h5")

#triqs_add_python_test(densdens)
#triqs_add_python_test(jperp)
#triqs_add_python_test(all)
#triqs_add_python_test(all_continuousBoson)

#set_property( DIRECTORY PROPERTY ENVIRONMENT PYTHONPATH=${CMAKE_BINARY_DIR}/python:$ENV{PYTHONPATH} )
add_test( NAME densdens COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/densdens.py )
set_property( TEST densdens PROPERTY ENVIRONMENT PYTHONPATH=${CMAKE_BINARY_DIR}/python:$ENV{PYTHONPATH} )
