# Build the python module
add_cpp2py_module(triqs_ctint)

target_link_libraries(triqs_ctint triqs_ctint_c)

# We need to include the convertes.hxx files
target_include_directories(triqs_ctint PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS triqs_ctint DESTINATION ${CPP2PY_PYTHON_LIB_DEST_ROOT})

# Add version file
configure_file(version.py.in version.py)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/version.py DESTINATION ${CPP2PY_PYTHON_LIB_DEST_ROOT})
