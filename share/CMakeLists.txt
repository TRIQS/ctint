add_subdirectory(cmake)

if(NOT CMAKE_INSTALL_PREFIX STREQUAL TRIQS_ROOT)

  configure_file(triqs_ctint.modulefile.in triqs_ctint.modulefile @ONLY)
  configure_file(triqs_ctintvars.sh.in triqs_ctintvars.sh @ONLY)

  install(
    FILES
      ${CMAKE_CURRENT_BINARY_DIR}/triqs_ctint.modulefile
      ${CMAKE_CURRENT_BINARY_DIR}/triqs_ctintvars.sh
    DESTINATION share
  )
 
  message(STATUS "***************************************************************")
  message(STATUS "* Custom install Location. Use: 			         ")
  message(STATUS "*                                                              ")
  message(STATUS "*   source ${CMAKE_INSTALL_PREFIX}/share/triqs_ctintvars.sh      ")
  message(STATUS "*                                                              ") 
  message(STATUS "* to set up the environment variables                          ") 
  message(STATUS "***************************************************************")

endif()
