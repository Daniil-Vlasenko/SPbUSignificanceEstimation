cmake_minimum_required(VERSION 3.22)
project(SPBUSignificanceEstimation)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(SPBUSignificanceEstimation SignificanceEstimation.cpp SignificanceEstimation.h
        probdist.cpp probdist.h main.cpp)
