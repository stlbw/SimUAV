cmake_minimum_required(VERSION 3.24)
project(SimUAV)

set(CMAKE_CXX_STANDARD 20)
add_subdirectory(Google_tests)
add_executable(SimUAV src/main.cpp declaredFun.h)
