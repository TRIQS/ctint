FILE(GLOB_RECURSE cpp_sources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)
add_library(ctint_c ${cpp_sources})

target_link_libraries(ctint_c PUBLIC triqs)
target_link_libraries(ctint_c PRIVATE nfft)

target_compile_options(ctint_c PUBLIC -std=c++17)

target_compile_definitions(ctint_c PUBLIC
			    CTINT_GIT_HASH=${CTINT_GIT_HASH}
			    TRIQS_GIT_HASH=${TRIQS_GIT_HASH}
			    $<$<CONFIG:Debug>:DEBUG_CTINT>
			    $<$<CONFIG:Debug>:TRIQS_ARRAYS_ENFORCE_BOUNDCHECK>
			    $<$<BOOL:${GTAU_IS_COMPLEX}>:GTAU_IS_COMPLEX>
			    $<$<BOOL:${INTERACTION_IS_COMPLEX}>:INTERACTION_IS_COMPLEX>
			  )

# Install
install(TARGETS ctint_c DESTINATION lib)
