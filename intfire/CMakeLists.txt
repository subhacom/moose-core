cmake_minimum_required(VERSION 3.0)
include( ${CMAKE_CURRENT_SOURCE_DIR}/../CheckCXXCompiler.cmake)
file(GLOB SRC *.cpp)
add_library(intfire ${SRC})
