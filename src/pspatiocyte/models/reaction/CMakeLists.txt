
set(method diffusion)

set(dest ${CMAKE_BINARY_DIR}/models/${method})
file(MAKE_DIRECTORY ${dest})

add_subdirectory(reversible)
#add_subdirectory(forward)


