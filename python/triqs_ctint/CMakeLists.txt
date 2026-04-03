# Configure the version
configure_file(version.py.in version.py)

# All Python files. Copy them in the build dir to have a complete package for the tests.
file(GLOB_RECURSE python_sources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.py)
file(GLOB_RECURSE wrap_generators RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *_desc.py)
list(REMOVE_ITEM python_sources "${wrap_generators}")
foreach(file ${python_sources})
  configure_file(${file} ${file} COPYONLY)
endforeach()

# Install python files to proper location
set(PYTHON_LIB_DEST ${TRIQS_PYTHON_LIB_DEST_ROOT}/${PROJECT_NAME})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/version.py DESTINATION ${PYTHON_LIB_DEST})
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} DESTINATION ${TRIQS_PYTHON_LIB_DEST_ROOT} FILES_MATCHING PATTERN "*.py" PATTERN "*_desc.py" EXCLUDE)

# Install triqs_hartree_fock dependency (pure Python)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/deps/triqs_hartree_fock/python/triqs_hartree_fock DESTINATION ${TRIQS_PYTHON_LIB_DEST_ROOT} FILES_MATCHING PATTERN "*.py")

# Build and install the real-variant python modules
foreach(gen ${wrap_generators})
  string(REPLACE "_desc.py" "" gen ${gen})
  get_filename_component(module_name ${gen} NAME_WE)
  get_filename_component(module_dir ${gen} DIRECTORY)
  add_cpp2py_module(NAME ${module_name} DIRECTORY ${module_dir})
  add_library(${PROJECT_NAME}::${module_name} ALIAS ${module_name})
  target_link_libraries(${module_name} ${PROJECT_NAME}_c triqs_py)
  install(TARGETS ${module_name} DESTINATION ${PYTHON_LIB_DEST}/${module_dir})
endforeach()

# ========= Complex Variant Modules ==========

# Function to build a variant of a cpp2py module from an already-generated _wrap.cxx.
# The variant library's PUBLIC compile definitions propagate automatically.
function(add_ctint_variant_module MODULE_NAME VARIANT_SUFFIX VARIANT_LIB)
  set(variant_target ${MODULE_NAME}_${VARIANT_SUFFIX})
  set(wrap_cxx ${CMAKE_CURRENT_BINARY_DIR}/${MODULE_NAME}_wrap.cxx)
  set(output_dir ${CMAKE_CURRENT_BINARY_DIR}/_${VARIANT_SUFFIX})

  add_library(${variant_target} MODULE ${wrap_cxx})
  target_link_libraries(${variant_target} ${VARIANT_LIB} triqs_py cpp2py::cpp2py)

  set_target_properties(${variant_target} PROPERTIES
    PREFIX ""
    OUTPUT_NAME ${MODULE_NAME}
    SUFFIX "${CPP2PY_PYTHON_MODULE_EXT}"
    LIBRARY_OUTPUT_DIRECTORY ${output_dir}
  )
  target_compile_options(${variant_target} PRIVATE
    $<$<PLATFORM_ID:Linux>:-fvisibility=hidden>
    $<$<PLATFORM_ID:Linux>:-fvisibility-inlines-hidden>
  )

  install(TARGETS ${variant_target} DESTINATION ${PYTHON_LIB_DEST}/_${VARIANT_SUFFIX})
endfunction()

# Helper to create a variant subpackage with __init__.py and all cpp2py modules
function(add_ctint_variant VARIANT_SUFFIX VARIANT_LIB)
  foreach(gen ${wrap_generators})
    string(REPLACE "_desc.py" "" gen ${gen})
    get_filename_component(module_name ${gen} NAME_WE)
    add_ctint_variant_module(${module_name} ${VARIANT_SUFFIX} ${VARIANT_LIB})
  endforeach()

  # Generate and install __init__.py for the subpackage
  set(init_file ${CMAKE_CURRENT_BINARY_DIR}/_${VARIANT_SUFFIX}/__init__.py)
  file(WRITE ${init_file} "")
  install(FILES ${init_file} DESTINATION ${PYTHON_LIB_DEST}/_${VARIANT_SUFFIX})
endfunction()

if(GTAU_IS_COMPLEX)
  add_ctint_variant(complex_gtau ${PROJECT_NAME}_c_complex_gtau)
endif()

if(INTERACTION_IS_COMPLEX)
  add_ctint_variant(complex_all ${PROJECT_NAME}_c_complex_all)
endif()
