set(method ode)
set(copy_files
  ode_distributive.csv
  ode_processive.csv
  ode_distributive.py
  ode_processive.py)

set(dest ${CMAKE_BINARY_DIR}/models/mapk/${method})
file(MAKE_DIRECTORY ${dest})

foreach(copy_file ${copy_files})
  get_filename_component(target ${copy_file} NAME_WE)
  add_custom_target(${target}_${method} ALL
    COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/${copy_file}
    ${dest})
endforeach(copy_file ${copy_file})


