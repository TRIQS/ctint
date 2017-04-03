# where will the python end up in triqs?
set(python_destination ${CMAKE_PROJECT_NAME})

# make a local pytriqs copy
triqs_prepare_local_pytriqs(${python_destination})

# Build the python module

triqs_python_extension_and_install(ctint ${python_destination})

target_link_libraries(ctint ctint_c)
include_directories(${TRIQS_INCLUDE_ALL} ${CMAKE_CURRENT_SOURCE_DIR})
triqs_set_rpath_for_target(ctint)


#triqs_python_extension_and_install(postproc)
#target_link_libraries(postproc ctint_c)
#include_directories(${TRIQS_INCLUDE_ALL} ${CMAKE_CURRENT_SOURCE_DIR})
#triqs_set_rpath_for_target(postproc)

# add version file

configure_file(version.py.in version.py)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/version.py DESTINATION ${TRIQS_PYTHON_LIB_DEST_ROOT}/${python_destination})
   

