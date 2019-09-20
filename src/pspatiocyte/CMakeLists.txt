
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(pspatiocyte_sources
  Compartment.cpp
  Lattice.cpp
  ParallelEnvironment.cpp
  Reaction.cpp
  Species.cpp
  World.cpp
  pSpatiocyte.cpp)

add_executable(pspatiocyte ${pspatiocyte_sources})

find_package(MPI REQUIRED)

find_package(Boost COMPONENTS random REQUIRED)

target_compile_options(pspatiocyte PRIVATE
  -O3 -mtune=native 
  ${MPI_CXX_COMPILE_FLAGS})

target_link_libraries(pspatiocyte PRIVATE 
  ${MPI_CXX_LIBRARIES}
  ${MPI_CXX_LINK_FLAGS}
  ${Boost_LIBRARIES})

target_include_directories(pspatiocyte PRIVATE
  ${MPI_CXX_INCLUDE_PATH}
  ${Boost_INCLUDE_DIRS})

