# specify CMake minimum version
cmake_minimum_required(VERSION 2.8)

# project name
project(EWlattice3)

# set the output folder where the program being created
# CMAKE_SOURCE_DIR: by default, it is the top-level source directory
# CMAKE_BINARY_DIR: top-level directory of the build tree
# EXECUTABLE_OUTPUT_PATH: by default, it is CMAKE_CURRENT_BINARY_DIR
# LIBRARY_OUTPUT_PATH: by default it is CMAKE_CURRENT_BINARY_DIR

set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# set CXX flags
set(CMAKE_CXX_FLAGS "-std=c++11 -m64 -O2 -ffast-math -flto -funroll-loops -march=native")

### external libraries

### MPI
find_package(MPI REQUIRED)
if (MPI_CXX_FOUND)
    include_directories(SYSTEM ${MPI_CXX_INCLUDE_DIRS})
    set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
    set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
    # COMPILER SETTING IN THIS WAY DOES NOT WORK
    # it must be set before any project().
    set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER}) 
else (MPI_CXX_FOUND)
    message(SEND_ERROR "This application cannot compile without MPI")
endif (MPI_CXX_FOUND)

### HDF5
find_package(HDF5 REQUIRED)
if (HDF5_FOUND)
    include_directories(SYSTEM ${HDF5_INCLUDE_DIRS})
else (HDF5_FOUND)
    message(SEND_ERROR "This application cannot compile without HDF5")
endif (HDF5_FOUND)
    
### FFTW
find_path(FFTW_INCLUDE_DIRS fftw3.h)
find_library(FFTW_LIBRARIES NAMES fftw3)
include_directories(SYSTEM ${FFTW_INCLUDE_DIRS})

### Eigen
# set CMAKE_MODULE_PATH and the find_package works in MODULE mode.
# the FindEigen3.cmake file is required in the CMAKE_MODULE_PATH folder.
find_package(Eigen3 REQUIRED MODULE)
if (EIGEN3_FOUND)
    include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})
else (EIGEN3_FOUND)
    message(SEND_ERROR "This application cannot compile without Eigen3")
endif (EIGEN3_FOUND)

# the following folders will be included
include_directories(${PROJECT_SOURCE_DIR})
include_directories(${ParaSite})
include_directories(${EW_Model})

# add source files
file(GLOB PARASITE_SRC "${PROJECT_SOURCE_DIR}/ParaSite/*.cpp")
file(GLOB EWMODEL_SRC "${PROJECT_SOURCE_DIR}/EW_Model/*.cpp")
set(SRCS "${PARASITE_SRC}" "${EWMODEL_SRC}" "${PROJECT_SOURCE_DIR}/main.cpp")

add_executable(emlat.out ${SRCS})
target_link_libraries(emlat.out ${MPI_CXX_LIBRARIES})
target_link_libraries(emlat.out ${HDF5_LIBRARIES})
target_link_libraries(emlat.out ${FFTW_LIBRARIES})



# print all variables
function(echo_all_cmake_variable_values)
  message(STATUS “”)
  get_cmake_property(vs VARIABLES)
  foreach(v ${vs})
    message(STATUS “${v}=’${${v}}'”)
  endforeach(v)
  message(STATUS “”)
endfunction()

# echo_all_cmake_variable_values()