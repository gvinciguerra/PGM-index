cmake_minimum_required(VERSION 3.10)

project(alglib)

FILE(GLOB alglib_HEADERS src/*.h)
FILE(GLOB alglib_SOURCES src/*.cpp)

add_library(alglib ${alglib_SOURCES})
target_include_directories(alglib PUBLIC src)