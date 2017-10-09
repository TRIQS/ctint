# Build the python module
add_cpp2py_module(ctint)

target_link_libraries(ctint ctint_c)

install(TARGETS ctint DESTINATION ${CPP2PY_PYTHON_LIB_DEST_ROOT}/triqs_ctint)

# Add version file
configure_file(version.py.in version.py)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/version.py DESTINATION ${CPP2PY_PYTHON_LIB_DEST_ROOT}/triqs_ctint)
