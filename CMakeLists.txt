cmake_minimum_required(VERSION 3.24)
project(SimUAV)

set(CMAKE_CXX_STANDARD 20)
add_subdirectory(Google_tests)
add_executable(SimUAV src/main.cpp declaredFun.h)

# Get the absolute path to the "database" directory
set(DATABASE_DIR "${CMAKE_SOURCE_DIR}/database")
set(PATH_DIR "${CMAKE_SOURCE_DIR}/path")
set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")
set(CMAKE_EXE_LINKER_FLAGS "-static")


# Add a compile-time definition for the database directory path
target_compile_definitions(SimUAV PRIVATE DATABASE_DIR="${DATABASE_DIR}")
target_compile_definitions(SimUAV PRIVATE PATH_DIR="${PATH_DIR}")
