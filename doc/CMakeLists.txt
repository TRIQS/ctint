# Generate the conf.py
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/conf.py.in ${CMAKE_CURRENT_BINARY_DIR}/conf.py @ONLY)

# ---------------------------------
# Top Sphinx target
# ---------------------------------
# Sources
file(GLOB_RECURSE sources *.rst)

set(sphinx_top ${CMAKE_CURRENT_BINARY_DIR}/html/contents.html)
add_custom_command(OUTPUT ${sphinx_top} DEPENDS ${sources}
  		   COMMAND ${TRIQS_SPHINXBUILD_EXECUTABLE} -c . -j8 -b html ${CMAKE_CURRENT_BINARY_DIR} html)
add_custom_target(docs_sphinx ALL DEPENDS ${sphinx_top} ${CMAKE_CURRENT_BINARY_DIR})
add_dependencies(docs_sphinx ctint) # py_copy

# ------------------------------------------------------------------------------------------------
# Create a temporary copy of the doc to be completed by the automatically generated one.
# ------------------------------------------------------------------------------------------------
#file(GLOB_RECURSE all_rst_files RELATIVE ${CMAKE_SOURCE_DIR}/doc *.rst _* *.py *.png *.cpp)
#set(rst_copy_tar  ${CMAKE_CURRENT_BINARY_DIR}/rst_copy.tar)
#add_custom_command (OUTPUT ${rst_copy_tar} DEPENDS ${all_rst_files}
		    #COMMAND cd ${CMAKE_SOURCE_DIR}/doc && tar cf ${rst_copy_tar} ${all_rst_files} && cd ${CMAKE_BINARY_DIR}/doc && tar xf ${rst_copy_tar} )
#add_custom_target(rst_copy DEPENDS ${rst_copy_tar})
#add_dependencies(docs_sphinx rst_copy) # we must first have copied the sources

# ---------------------------------
# Run c++2rst
# ---------------------------------
add_custom_command (OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/cpp2rst.log
                    COMMAND c++2rst -N triqs_ctint
		    ${CMAKE_SOURCE_DIR}/c++/solver_core.hpp 
		    --output_directory=${CMAKE_CURRENT_BINARY_DIR}/cpp2rst_generated
		    --includes=${CMAKE_SOURCE_DIR}/c++ 
		    --includes=${TRIQS_ROOT}/include 
		    2>&1 > cpp2rst.log
		    ) 
add_custom_target(docs_cpp2rst DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/cpp2rst.log) 
add_dependencies(docs_cpp2rst ctint)
add_dependencies(docs_sphinx docs_cpp2rst)

# ---------------------------------
# Install
# ---------------------------------
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html/ COMPONENT documentation DESTINATION share/doc/ctint
 FILES_MATCHING 
 PATTERN "*.html"  # Use REGEX ??
 PATTERN "*.png"
 PATTERN "*.jpg"
 PATTERN "*.js"
 PATTERN "_*"
 PATTERN "*.png" 
 PATTERN "*.gif" 
 PATTERN "*.xsl" 
 PATTERN "*.css"
 PATTERN "*.pdf"
 PATTERN "*.py"
 PATTERN "*.txt"
 PATTERN "*.inv"
 PATTERN "*.bib"
 )

