cmake_minimum_required(VERSION 3.10)
project(REMEZ_REBOOT)

add_definitions("-g -Wall -march=native -lntl -lgmp -lm")
add_executable(a.out main.cpp Point.cpp func.cpp Polynomial.cpp RemezGeneral.cpp)

target_compile_options(a.out PUBLIC -Wall)
target_include_directories(a.out PUBLIC ${CMAKE_SOURCE_DIR}/includes)

target_link_libraries(a.out ntl gmp m)
