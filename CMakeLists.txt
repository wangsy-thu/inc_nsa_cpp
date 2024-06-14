cmake_minimum_required(VERSION 3.28)
project(inc_nas_cpp)

set(CMAKE_CXX_STANDARD 17)

add_executable(inc_nas_cpp main.cpp
        graph/Graph.cpp
        graph/Graph.h
        community/DynComDetector.cpp
        community/DynComDetector.h)
