cmake_minimum_required(VERSION 3.19)
project(bipartite_match)

set(CMAKE_CXX_STANDARD 17)

add_executable(bipartite_match main.cpp Hungarian.cpp bipartite_match.cpp bipartite_match.h)