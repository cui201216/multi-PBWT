cmake_minimum_required(VERSION 3.28)
project(multiPBWT)

set(CMAKE_CXX_STANDARD  14)

# 设置静态链接标志
set(CMAKE_EXE_LINKER_FLAGS "-static")

add_executable(multiPBWT main.cpp)
