cmake_minimum_required(VERSION 2.8.0 FATAL_ERROR)

set(CMAKE_CXX_COMPILER "icpc" CACHE STRING "" FORCE)

set(GTEST_ROOT $ENV{HOME}/usr/src/gtest-1.7.0)
include_directories(
    ${GTEST_ROOT}/include
    )
add_executable(gtest_example gtest_example.cpp)
target_link_libraries(gtest_example
    pthread
    ${GTEST_ROOT}/build/libgtest.a
    ${GTEST_ROOT}/build/libgtest_main.a
    )

