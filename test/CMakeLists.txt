# include guard
cmake_minimum_required(VERSION 2.8.0 )

add_test(
  NAME lanczos_hubbard_square
  COMMAND ${CMAKE_SOURCE_DIR}/test/lanczos_hubbard_square.sh
)
