# generate the conf.py
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/conf.py.in ${CMAKE_CURRENT_BINARY_DIR}/conf.py)

# create documentation target
set(DOC_SOURCE ${CMAKE_CURRENT_BINARY_DIR})
set(sphinx_top ${CMAKE_CURRENT_BINARY_DIR}/html/index.html)
add_custom_command(OUTPUT ${sphinx_top} DEPENDS py_copy ctint ctint_c
                   COMMAND ${CMAKE_BINARY_DIR}/build_pytriqs ${TRIQS_SPHINXBUILD_EXECUTABLE} -c . -b html ${DOC_SOURCE} html)
add_custom_target(doc_sphinx ALL DEPENDS ${sphinx_top} ${DOC_SOURCE})

# Create a temporary copy of the doc that will be completed by the automatically generated one.
FILE(GLOB_RECURSE all_rst_files RELATIVE ${CMAKE_SOURCE_DIR}/doc *.rst _* *.py *.png *.cpp)
SET(rst_copy_tar  ${CMAKE_CURRENT_BINARY_DIR}/rst_copy.tar)
add_custom_command (OUTPUT ${rst_copy_tar} DEPENDS ${all_rst_files}
 COMMAND cd ${CMAKE_SOURCE_DIR}/doc && tar cf ${rst_copy_tar} ${all_rst_files} && cd ${CMAKE_BINARY_DIR}/doc && tar xf ${rst_copy_tar} )
add_custom_target(rst_copy DEPENDS ${rst_copy_tar})
add_dependencies(doc_sphinx rst_copy) # we must first have copied the sources

# use c++2doc to automatically build the documentation of C++  
if (DocWithCpp2doc)
 FILE(GLOB_RECURSE all_header_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ../c++/*.hpp)
 add_custom_command (OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/cpp2py.log DEPENDS ${all_header_files}  COMMAND export PYTHONPATH=\$PYTHONPATH:${CMAKE_INSTALL_PREFIX}/share/triqs/cpp2doc/mako && ${PYTHON_INTERPRETER} ${CMAKE_INSTALL_PREFIX}/bin/c++2doc.py -N triqs -N applications -N impurity_solvers ${CMAKE_CURRENT_SOURCE_DIR}/common.hpp 
      --output_directory=${CMAKE_CURRENT_BINARY_DIR}/reference
      --includes=${CMAKE_SOURCE_DIR}/c++ --includes=${CMAKE_INSTALL_PREFIX}/include --includes=${CMAKE_INSTALL_PREFIX}/share/triqs/cmake  
      --mako_location=${CMAKE_INSTALL_PREFIX}/share/triqs 2>&1 > cpp2py.log) 
add_custom_target(docs_cpp2doc DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/cpp2py.log) 
add_dependencies(doc_sphinx docs_cpp2doc) # we must have installed the sources before building
endif()

# install
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html/ COMPONENT documentation DESTINATION share/doc/ctint
 FILES_MATCHING 
 PATTERN "*.html" 
 PATTERN "*.png"
 PATTERN "*.js"
 PATTERN "_*"
 PATTERN "*.png" 
 PATTERN "*.jpg" 
 PATTERN "*.gif" 
 PATTERN "*.xsl" 
 PATTERN "*.css"
 PATTERN "*.pdf"
 PATTERN "*.py"
 PATTERN "*.txt"
 PATTERN "*.inv"
 PATTERN "*.bib"
 )

