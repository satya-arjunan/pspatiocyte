
set(method readdy_serial)
set(copy_files
  model.py
  run.py
  elapsed_time.txt
  output.txt
  plot_output.py)

set(dest ${CMAKE_BINARY_DIR}/models/benchmark/${method})
file(MAKE_DIRECTORY ${dest})

foreach(copy_file ${copy_files})
  get_filename_component(target ${copy_file} NAME_WE)
  add_custom_target(${target}_${method} ALL
    COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/${copy_file}
    ${dest})
endforeach(copy_file ${copy_file})


