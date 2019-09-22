
set(method benchmark)
set(copy_files
  plot_benchmark.py
  plot_output.py)

set(dest ${CMAKE_BINARY_DIR}/models/${method})
file(MAKE_DIRECTORY ${dest})

foreach(copy_file ${copy_files})
  get_filename_component(target ${copy_file} NAME_WE)
  add_custom_target(${target}_${method} ALL
    COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/${copy_file}
    ${dest})
endforeach(copy_file ${copy_file})

add_subdirectory(pspatiocyte)
add_subdirectory(spatiocyte)
add_subdirectory(smoldyn)
add_subdirectory(readdy)
add_subdirectory(readdy_serial)


