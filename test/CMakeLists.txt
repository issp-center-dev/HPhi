# include guard
cmake_minimum_required(VERSION 2.8.0 )

#
# Lanczos
#
add_test(
  NAME lanczos_hubbard_square
  COMMAND ${CMAKE_SOURCE_DIR}/test/lanczos_hubbard_square.sh
)
add_test(
  NAME lanczos_kondo_chain
  COMMAND ${CMAKE_SOURCE_DIR}/test/lanczos_kondo_chain.sh
)
add_test(
  NAME lanczos_spin_kagome
  COMMAND ${CMAKE_SOURCE_DIR}/test/lanczos_spin_kagome.sh
)
add_test(
  NAME lanczos_genspin_ladder
  COMMAND ${CMAKE_SOURCE_DIR}/test/lanczos_genspin_ladder.sh
)
add_test(
  NAME lanczos_hubbardgc_tri
  COMMAND ${CMAKE_SOURCE_DIR}/test/lanczos_hubbardgc_tri.sh
)
add_test(
  NAME lanczos_kondogc_chain
  COMMAND ${CMAKE_SOURCE_DIR}/test/lanczos_kondogc_chain.sh
)
add_test(
  NAME lanczos_spingc_honey
  COMMAND ${CMAKE_SOURCE_DIR}/test/lanczos_spingc_honey.sh
)
add_test(
  NAME lanczos_genspingc_ladder
  COMMAND ${CMAKE_SOURCE_DIR}/test/lanczos_genspingc_ladder.sh
)
#
# LOBCG
#
add_test(
  NAME lobcg_hubbard_square
  COMMAND ${CMAKE_SOURCE_DIR}/test/lobcg_hubbard_square.sh
)
add_test(
  NAME lobcg_kondo_chain
  COMMAND ${CMAKE_SOURCE_DIR}/test/lobcg_kondo_chain.sh
)
add_test(
  NAME lobcg_spin_kagome
  COMMAND ${CMAKE_SOURCE_DIR}/test/lobcg_spin_kagome.sh
)
add_test(
  NAME lobcg_genspin_ladder
  COMMAND ${CMAKE_SOURCE_DIR}/test/lobcg_genspin_ladder.sh
)
add_test(
  NAME lobcg_hubbardgc_tri
  COMMAND ${CMAKE_SOURCE_DIR}/test/lobcg_hubbardgc_tri.sh
)
add_test(
  NAME lobcg_kondogc_chain
  COMMAND ${CMAKE_SOURCE_DIR}/test/lobcg_kondogc_chain.sh
)
add_test(
  NAME lobcg_spingc_honey
  COMMAND ${CMAKE_SOURCE_DIR}/test/lobcg_spingc_honey.sh
)
add_test(
  NAME lobcg_genspingc_ladder
  COMMAND ${CMAKE_SOURCE_DIR}/test/lobcg_genspingc_ladder.sh
)
#
# Full Diag
#
add_test(
  NAME fulldiag_hubbard_chain
  COMMAND ${CMAKE_SOURCE_DIR}/test/fulldiag_hubbard_chain.sh
)
add_test(
  NAME fulldiag_kondo_chain
  COMMAND ${CMAKE_SOURCE_DIR}/test/fulldiag_kondo_chain.sh
)
add_test(
  NAME fulldiag_spin_tri
  COMMAND ${CMAKE_SOURCE_DIR}/test/fulldiag_spin_tri.sh
)
add_test(
  NAME fulldiag_genspin_ladder
  COMMAND ${CMAKE_SOURCE_DIR}/test/fulldiag_genspin_ladder.sh
)
add_test(
  NAME fulldiag_hubbardgc_tri
  COMMAND ${CMAKE_SOURCE_DIR}/test/fulldiag_hubbardgc_tri.sh
)
add_test(
  NAME fulldiag_kondogc_chain
  COMMAND ${CMAKE_SOURCE_DIR}/test/fulldiag_kondogc_chain.sh
)
add_test(
  NAME fulldiag_spingc_tri
  COMMAND ${CMAKE_SOURCE_DIR}/test/fulldiag_spingc_tri.sh
)
add_test(
  NAME fulldiag_genspingc_ladder
  COMMAND ${CMAKE_SOURCE_DIR}/test/fulldiag_genspingc_ladder.sh
)
