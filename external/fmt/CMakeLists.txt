cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

add_library(fmt STATIC src/format.cc src/os.cc)
target_include_directories(fmt PRIVATE include)
