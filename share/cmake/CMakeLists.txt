configure_file(triqs_ctint-config.cmake.in triqs_ctint-config.cmake @ONLY)
configure_file(triqs_ctint-config-version.cmake.in triqs_ctint-config-version.cmake @ONLY)
install(
  FILES
    ${CMAKE_CURRENT_BINARY_DIR}/triqs_ctint-config.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/triqs_ctint-config-version.cmake
  DESTINATION lib/cmake/triqs_ctint
)

install(EXPORT ctint-targets DESTINATION lib/cmake/triqs_ctint)
