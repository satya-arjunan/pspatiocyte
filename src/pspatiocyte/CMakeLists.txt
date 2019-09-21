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
find_package(Boost COMPONENTS random REQUIRED)

target_compile_options(${PROJECT_NAME} PRIVATE
  -O3 -mtune=native 
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

file(GLOB MODEL_DIRS models/*)
foreach(model_dir ${MODEL_DIRS})
  get_filename_component(model_name ${model_dir} NAME_WE)
  set(model_src ${model_dir}/${model_name}.cpp)
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
endforeach(model_dir ${MODEL_DIRS})
