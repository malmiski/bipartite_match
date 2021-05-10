cmake_minimum_required(VERSION 3.19)
project(bipartite_match)

set(CMAKE_CXX_STANDARD 17)
add_executable(bipartite_match main.cpp Hungarian.cpp bipartite_match.cpp bipartite_match.h  wspd.cpp wspd.h bipartite_match_ANN.cpp bipartite_match_ANN.h)
target_link_libraries(bipartite_match ${CMAKE_SOURCE_DIR}/lib/libANN.a)
