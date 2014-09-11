cmake_minimum_required(VERSION 2.8.4)
project(DSEpp)

#Set CFLAGs
ADD_DEFINITIONS(
    -fopenmp
    -O3
    -std=gnu++0x
)

#Set LFLAGs
set (CMAKE_EXE_LINKER_FLAGS -fopenmp )

#define path to external libraries
include_directories(Libraries/qft++/include)
include_directories(Libraries/eigen)

#grab all file we need to compile
file(GLOB_RECURSE SRC
    "source/*/*.cpp"
)

add_library( myLib ${SRC} )
target_include_directories (myLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

set(SOURCE_FILES source/main.cpp)

add_executable(DSEpp ${SOURCE_FILES})
target_link_libraries( DSEpp myLib )