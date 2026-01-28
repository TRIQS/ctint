file(GLOB_RECURSE all_benchs RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

foreach(bench ${all_benchs})
  get_filename_component(bench_name ${bench} NAME_WE)
  add_executable(${bench_name} ${bench})
  target_link_libraries(${bench_name} ${PROJECT_NAME}_c finufft benchmark_main)
endforeach()
