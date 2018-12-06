# Build the python modules
add_cpp2py_module(solver_core)
add_cpp2py_module(post_process)

target_link_libraries(solver_core ctint_c)
target_link_libraries(post_process ctint_c)

# We need to include the convertes.hxx files
target_include_directories(solver_core PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(post_process PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# Configure the version
configure_file(version.py.in version.py)

# All Python files. Copy them in the build dir to have a complete package for the tests.
set(PYTHON_SOURCES __init__.py solver.py)
foreach(f ${PYTHON_SOURCES})
  configure_file(${f} ${f} COPYONLY)
endforeach()

# Install python module to proper location
set(PYTHON_LIB_DEST ${CPP2PY_PYTHON_LIB_DEST_ROOT}/triqs_ctint)
install(TARGETS solver_core post_process DESTINATION ${PYTHON_LIB_DEST})
install(FILES ${PYTHON_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/version.py DESTINATION ${PYTHON_LIB_DEST})
