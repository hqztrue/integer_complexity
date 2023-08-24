# this file is used by some software that help the developers.
# you can just ignore this file if you feel that's better.

# manual DEVELOPMENT compilation is done with "gcc -Wall -pedantic -O0 main.c -o qs"
# manual PRODUCTION compilation is done with "gcc -Wall -pedantic -O3 main.c -o qs"

cmake_minimum_required(VERSION 3.20)
project(qs C)
set(CMAKE_C_STANDARD 99)
set (CMAKE_C_FLAGS "-Wall -Wextra -pedantic -g -O3 -std=c99")
add_executable(qs main.c)
