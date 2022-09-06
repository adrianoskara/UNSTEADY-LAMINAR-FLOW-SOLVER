cmake_minimum_required(VERSION 3.22)
project(unsteady)

set(CMAKE_CXX_STANDARD 14)
include_directories('/usr/include/eigen-3.4.0/')
include_directories('/home/tkaravasilhs/CLionProjects/unsteady/progressbar-master')

add_executable(unsteady main.cpp Variables.h Variables.cpp schemes.cpp schemes.h)
