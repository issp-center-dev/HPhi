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
  NAME lanczos_spin_ladder_DM
COMMAND ${CMAKE_SOURCE_DIR}/test/lanczos_spin_ladder_DM.sh  
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
  NAME lanczos_spingc_hcor
  COMMAND ${CMAKE_SOURCE_DIR}/test/lanczos_spingc_hcor.sh
)
add_test(
  NAME lanczos_genspingc_ladder
  COMMAND ${CMAKE_SOURCE_DIR}/test/lanczos_genspingc_ladder.sh
)
add_test(
  NAME lanczos_hubbard_square_restart
  COMMAND ${CMAKE_SOURCE_DIR}/test/lanczos_hubbard_square_restart.sh
)
add_test(
  NAME lanczos_kondo_chain_restart
  COMMAND ${CMAKE_SOURCE_DIR}/test/lanczos_kondo_chain_restart.sh
)
add_test(
  NAME lanczos_spin_kagome_restart
  COMMAND ${CMAKE_SOURCE_DIR}/test/lanczos_spin_kagome_restart.sh
)
add_test(
  NAME lanczos_genspin_ladder_restart
  COMMAND ${CMAKE_SOURCE_DIR}/test/lanczos_genspin_ladder_restart.sh
)
#
# LOBCG
#
add_test(
  NAME lobcg_hubbard_square
  COMMAND ${CMAKE_SOURCE_DIR}/test/lobcg_hubbard_square.sh
)
add_test(
  NAME lobcg_hubbard_szunconserv
  COMMAND ${CMAKE_SOURCE_DIR}/test/lobcg_hubbard_szunconserv.sh
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
  NAME lobcg_spingc_Sy
  COMMAND ${CMAKE_SOURCE_DIR}/test/lobcg_spingc_Sy.sh
)
add_test(
  NAME lobcg_spingc_hcor
  COMMAND ${CMAKE_SOURCE_DIR}/test/lobcg_spingc_hcor.sh
)
add_test(
  NAME lobcg_spingc_6body
  COMMAND ${CMAKE_SOURCE_DIR}/test/lobcg_spingc_6body.sh
)
add_test(
  NAME lobcg_genspingc_ladder
  COMMAND ${CMAKE_SOURCE_DIR}/test/lobcg_genspingc_ladder.sh
)
add_test(
  NAME lobcg_hubbard_square_restart
  COMMAND ${CMAKE_SOURCE_DIR}/test/lobcg_hubbard_square_restart.sh
)
add_test(
  NAME lobcg_kondo_chain_restart
  COMMAND ${CMAKE_SOURCE_DIR}/test/lobcg_kondo_chain_restart.sh
)
add_test(
  NAME lobcg_spin_kagome_restart
  COMMAND ${CMAKE_SOURCE_DIR}/test/lobcg_spin_kagome_restart.sh
)
add_test(
  NAME lobcg_genspin_ladder_restart
  COMMAND ${CMAKE_SOURCE_DIR}/test/lobcg_genspin_ladder_restart.sh
)
#
# Full Diag
#
add_test(
  NAME fulldiag_hubbard_chain
  COMMAND ${CMAKE_SOURCE_DIR}/test/fulldiag_hubbard_chain.sh
)
add_test(
  NAME fulldiag_Hamiltonian_IO
  COMMAND ${CMAKE_SOURCE_DIR}/test/fulldiag_ham_io.sh
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
#
# Spectrum
#
add_test(
  NAME spectrum_hubbard_square
  COMMAND ${CMAKE_SOURCE_DIR}/test/spectrum_hubbard_square.sh
)
add_test(
  NAME spectrum_kondo_chain
  COMMAND ${CMAKE_SOURCE_DIR}/test/spectrum_kondo_chain.sh
)
add_test(
  NAME spectrum_spin_kagome
  COMMAND ${CMAKE_SOURCE_DIR}/test/spectrum_spin_kagome.sh
)
add_test(
  NAME spectrum_genspin_ladder
  COMMAND ${CMAKE_SOURCE_DIR}/test/spectrum_genspin_ladder.sh
)
add_test(
  NAME spectrum_hubbardgc_tri
  COMMAND ${CMAKE_SOURCE_DIR}/test/spectrum_hubbardgc_tri.sh
)
add_test(
  NAME spectrum_kondogc_chain
  COMMAND ${CMAKE_SOURCE_DIR}/test/spectrum_kondogc_chain.sh
)
add_test(
  NAME spectrum_spingc_honey
  COMMAND ${CMAKE_SOURCE_DIR}/test/spectrum_spingc_honey.sh
)

#
# Time evolution
#
add_test(
  NAME te_ac_hubbard_square
  COMMAND ${CMAKE_SOURCE_DIR}/test/te_ac_hubbard_square.sh
)
add_test(
  NAME te_dc_hubbard_square
  COMMAND ${CMAKE_SOURCE_DIR}/test/te_dc_hubbard_square.sh
)
add_test(
  NAME te_pulse_hubbard_square
  COMMAND ${CMAKE_SOURCE_DIR}/test/te_pulse_hubbard_square.sh
)
add_test(
  NAME te_quench_hubbard_square
  COMMAND ${CMAKE_SOURCE_DIR}/test/te_quench_hubbard_square.sh
)
add_test(
  NAME te_hubbard_chain_interall
  COMMAND ${CMAKE_SOURCE_DIR}/test/te_hubbard_chain_interall.sh ${CMAKE_SOURCE_DIR}
)
add_test(
  NAME te_hubbard_chain_interall_diagonal
  COMMAND ${CMAKE_SOURCE_DIR}/test/te_hubbard_chain_interall_diagonal.sh ${CMAKE_SOURCE_DIR}
)
add_test(
  NAME te_spin_chain_interall
  COMMAND ${CMAKE_SOURCE_DIR}/test/te_spin_chain_interall.sh ${CMAKE_SOURCE_DIR}
)
add_test(
  NAME te_kondo_chain_interall
  COMMAND ${CMAKE_SOURCE_DIR}/test/te_kondo_chain_interall.sh ${CMAKE_SOURCE_DIR}
)
#
#mTPQ
#
add_test(
  NAME tpq_spin_kagome
  COMMAND ${CMAKE_SOURCE_DIR}/test/tpq_spin_kagome.sh ${CMAKE_SOURCE_DIR}
)
add_test(
  NAME tpq_spin_kagome_randomsphere
  COMMAND ${CMAKE_SOURCE_DIR}/test/tpq_spin_kagome_randomsphere.sh ${CMAKE_SOURCE_DIR}
)

#
#cTPQ
#
add_test(
  NAME ctpq_spin_kagome
  COMMAND ${CMAKE_SOURCE_DIR}/test/ctpq_spin_kagome.sh ${CMAKE_SOURCE_DIR}
)
add_test(
  NAME ctpq_spin_kagome_randomsphere
  COMMAND ${CMAKE_SOURCE_DIR}/test/ctpq_spin_kagome_randomsphere.sh ${CMAKE_SOURCE_DIR}
)


