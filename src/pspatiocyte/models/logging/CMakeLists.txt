
set(model_name logging)

set(copy_files
  run.py
  log_hostfile)

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
  add_custom_command(
    TARGET ${model_name} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy 
    ${CMAKE_CURRENT_SOURCE_DIR}/${copy_file}
    $<TARGET_FILE_DIR:${model_name}>)
endforeach()

