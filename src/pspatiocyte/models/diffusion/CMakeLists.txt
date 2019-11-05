
set(method diffusion)

set(dest ${CMAKE_BINARY_DIR}/models/${method})
file(MAKE_DIRECTORY ${dest})

add_subdirectory(trajectory)
add_subdirectory(dilute)
add_subdirectory(crowded)


