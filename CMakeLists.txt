# specify CMake minimum version
cmake_minimum_required(VERSION 2.8)

# project name
project(EWlattice3)
set(TARGET emlat.out)

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
#if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
#	set(CMAKE_CXX_FLAGS "-std=c++11 -m64 -O3 -ffast-math -flto -funroll-loops -march=native -Wall")
#elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
#	set(CMAKE_CXX_FLAGS "-std=c++11 -O3 -no-prec-div -ipo -xcore-avx2")
#endif()
set(CMAKE_CXX_FLAGS "-std=c++11 -m64 -O3 -ffast-math -flto -funroll-loops -march=native -Wall")

### external libraries

### MPI
find_package(MPI REQUIRED)
if (MPI_CXX_FOUND)
    include_directories(SYSTEM ${MPI_CXX_INCLUDE_DIRS})
    set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
    set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
    # COMPILER SETTING IN THIS WAY DOES NOT WORK
    # it must be set before any project().

    # if use Intel compiler, then use mpiicpc
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    	set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
	    get_filename_component(MPI_CXX_COMPILER_DIR ${MPI_CXX_COMPILER} DIRECTORY)
	    set(CMAKE_CXX_COMPILER "${MPI_CXX_COMPILER_DIR}/mpiicpc")
    endif()
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

add_executable(${TARGET} ${SRCS})
target_link_libraries(${TARGET} ${MPI_CXX_LIBRARIES})
target_link_libraries(${TARGET} ${HDF5_LIBRARIES})
target_link_libraries(${TARGET} ${FFTW_LIBRARIES})

### Intel MKL (optional)
find_package(MKL MODULE)
if (MKL_FOUND)
    include_directories(${MKL_INCLUDE_DIRS})
    target_link_libraries(${TARGET} ${MKL_LIBRARIES})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DEIGEN_USE_MKL_ALL")
endif()


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