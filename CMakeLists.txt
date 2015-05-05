cmake_minimum_required(VERSION 2.8.12)
project(DSEpp)

#Set CFLAGs
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fopenmp -g -std=gnu++0x")

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fopenmp -std=gnu++0x -ffast-math -O3")


#Set LFLAGs
set (CMAKE_EXE_LINKER_FLAGS -fopenmp )

#define path to external libraries
include_directories(Libraries/qft++/include)
include_directories(Libraries/eigen)

#grab all files we need to compile
file(GLOB_RECURSE SRC
    "source/*/*.cpp"
    #"source/*/*.hpp"
    #"source/*/*.h"
    #"source/*.hpp"
    #"source/*.h"
)

add_library( myLib ${SRC} )
target_include_directories (myLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

set(SOURCE_FILES source/Example.cpp)

add_executable(DSEpp ${SOURCE_FILES})
target_link_libraries( DSEpp myLib )