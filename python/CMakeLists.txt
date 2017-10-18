# Build the python module
add_cpp2py_module(solver_core)

target_link_libraries(solver_core triqs_ctint_c)

# We need to include the convertes.hxx files
target_include_directories(solver_core PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS solver_core DESTINATION ${CPP2PY_PYTHON_LIB_DEST_ROOT}/triqs_ctint)

# Add version file
configure_file(version.py.in version.py)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/version.py DESTINATION ${CPP2PY_PYTHON_LIB_DEST_ROOT}/triqs_ctint)
