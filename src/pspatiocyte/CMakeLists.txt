cmake_minimum_required(VERSION 3.9)
project(pSpatiocyte VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(GNUInstallDirs)

add_library(${PROJECT_NAME} SHARED
  src/Compartment.cpp
  src/Lattice.cpp
  src/ParallelEnvironment.cpp
  src/Reaction.cpp
  src/Species.cpp
  src/World.cpp)

find_package(MPI REQUIRED)
find_package(Boost COMPONENTS random filesystem REQUIRED)

target_compile_options(${PROJECT_NAME} PRIVATE
  #-O3 -mtune=native 
  -O3 -g3 -mtune=native # for valgrind debug
  ${MPI_CXX_COMPILE_FLAGS})

target_include_directories(${PROJECT_NAME} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    PRIVATE src
    ${MPI_CXX_INCLUDE_PATH}
    ${Boost_INCLUDE_DIRS})

target_link_libraries(${PROJECT_NAME} PRIVATE 
  ${MPI_CXX_LIBRARIES}
  ${MPI_CXX_LINK_FLAGS}
  ${Boost_LIBRARIES})

set_target_properties(${PROJECT_NAME} PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION 1)

install(TARGETS ${PROJECT_NAME} EXPORT ${PROJECT_NAME}Config
    ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR})
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME})

install(EXPORT ${PROJECT_NAME}Config DESTINATION share/${PROJECT_NAME}/cmake)

export(TARGETS ${PROJECT_NAME} FILE ${PROJECT_NAME}Config.cmake)

set(script_files
  gather_timecourse.py
  gather_coordinates.py)

foreach(script_file ${script_files})
  add_custom_command(
    TARGET ${PROJECT_NAME} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy 
    ${CMAKE_CURRENT_SOURCE_DIR}/scripts/${script_file}
    ${CMAKE_BINARY_DIR}/scripts/${script_file})
endforeach()

add_subdirectory(models/benchmark)
add_subdirectory(models/mapk)
add_subdirectory(models/trajectory)
add_subdirectory(models/dissociation)
add_subdirectory(models/diffusion)
add_subdirectory(models/reaction)
