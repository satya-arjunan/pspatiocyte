
set(method reaction)

set(dest ${CMAKE_BINARY_DIR}/models/${method})
file(MAKE_DIRECTORY ${dest})

add_subdirectory(reversible)
add_subdirectory(irreversible)
add_subdirectory(crowded_dissociation)


