cmake_minimum_required(VERSION 3.12)
project(PiecewiseGeometricModelIndex
        VERSION 1.0
        HOMEPAGE_URL https://github.com/gvinciguerra/PGM-index)

option(BUILD_EXAMPLES "Build the examples" ON)
option(BUILD_PGM_TUNER "Build the tuner target" ON)
option(BUILD_PGM_BENCHMARK "Build the benchmark target" ON)

set(CMAKE_CXX_STANDARD 17)

include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
if (COMPILER_SUPPORTS_MARCH_NATIVE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif ()

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
endif ()

find_package(OpenMP)
if (OpenMP_CXX_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    link_libraries(OpenMP::OpenMP_CXX)
endif ()

include_directories(include)

# PGM-index library
add_library(pgmindexlib INTERFACE)
target_include_directories(pgmindexlib INTERFACE include/)

if (BUILD_PGM_TUNER)
    add_subdirectory(tuner)
endif ()

if (BUILD_EXAMPLES)
    add_subdirectory(examples)
endif ()

if (BUILD_PGM_BENCHMARK)
    add_subdirectory(benchmark)
endif ()

add_subdirectory(c-interface)

enable_testing()
add_subdirectory(test)