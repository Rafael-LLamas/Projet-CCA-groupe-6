cmake_minimum_required(VERSION 3.10.0)
project(model_project VERSION 0.1.0 LANGUAGES C)

include(CTest)
enable_testing()

add_compile_options(-g -O0 -Wall)
set(CMAKE_C_STANDARD 99)

