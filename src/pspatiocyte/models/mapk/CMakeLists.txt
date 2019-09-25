
set(model_name mapk)

set(copy_files
  plot_ratio.py
  plot_output.py)

set(model_src ${model_name}.cpp)
add_executable(${model_name} ${model_src})
target_link_libraries(${model_name} PRIVATE 
  ${PROJECT_NAME}
  ${MPI_CXX_LIBRARIES}
  ${MPI_CXX_LINK_FLAGS}
  ${Boost_LIBRARIES})
target_include_directories(${model_name} PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  ${MPI_CXX_INCLUDE_PATH}
  ${Boost_INCLUDE_DIRS})
set_target_properties(${model_name}
  PROPERTIES
  ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
  LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/models/${model_name}")

foreach(copy_file ${copy_files})
  get_filename_component(target ${copy_file} NAME_WE)
  add_custom_target(${target}_${method} ALL
    COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/${copy_file}
    ${dest})
endforeach(copy_file ${copy_file})
