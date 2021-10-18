### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# CMake configuration
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
cmake_minimum_required(VERSION 3.16)

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Project configuration
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
project(ADTWNN1 VERSION 0.0.1)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
set(CMAKE_VERBOSE_MAKEFILE ON)

# Debug build: Activate GLIBC assertions
add_compile_definitions($<$<CONFIG:DEBUG>:_GLIBCXX_ASSERTIONS>)

find_package(Threads REQUIRED)

add_executable(ADTWNN1
        src/utils.hpp
        src/main.cpp
        src/tseries/readers/tsreader/tsreader.hpp
        src/tseries/readers/tsreader/tsreader.cpp
        src/tseries/tseries.hpp
    )

target_link_libraries(ADTWNN1 PRIVATE Threads::Threads)
