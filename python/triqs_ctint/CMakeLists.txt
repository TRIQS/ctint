# Build the python module
add_cpp2py_module(solver_core)

target_link_libraries(solver_core triqs_ctint_c)

# We need to include the convertes.hxx files
target_include_directories(solver_core PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# Add version file
configure_file(version.py.in version.py)

set(PYTHON_SOURCES
 ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py
 ${CMAKE_CURRENT_BINARY_DIR}/version.py
 )

# Install python module to proper location
set(PYTHON_LIB_DEST ${CPP2PY_PYTHON_LIB_DEST_ROOT}/triqs_ctint)
install(TARGETS solver_core DESTINATION ${PYTHON_LIB_DEST})
install(FILES ${PYTHON_SOURCES} DESTINATION ${PYTHON_LIB_DEST})
