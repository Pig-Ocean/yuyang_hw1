cmake_minimum_required(VERSION 3.10)
project(hw1 VERSION 1.0)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable(hw1
    src/main.c
    src/algebra.c)
target_include_directories(hw1 PRIVATE inc)