cmake_minimum_required(VERSION 3.20)
project(sample)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

add_executable(sample geomutils.cpp convhull.cpp)