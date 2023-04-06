cmake_minimum_required(VERSION 3.24)
project(SimUAV)

set(CMAKE_CXX_STANDARD 20)

add_executable(SimUAV main.cpp declaredFun.h Separated_reading_function/read_engine.cpp Separated_reading_function/read_propeller.cpp preprocessing/readingBatEnProp.cpp)
