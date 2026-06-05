# include guard for CMakeLists.txt
if(${CMAKE_PROJECT_NAME} STREQUAL "Project")
    message(FATAL_ERROR "cmake should not be executed directly in subdirectory.")
endif(${CMAKE_PROJECT_NAME} STREQUAL "Project")

# Macro for adding standard tests
macro(add_hphi_test test_name)
    add_test(
        NAME ${test_name}
        COMMAND ${CMAKE_SOURCE_DIR}/test/${test_name}.sh
    )
endmacro()

# Macro for adding tests with CMAKE_SOURCE_DIR argument
macro(add_hphi_test_with_srcdir test_name)
    add_test(
        NAME ${test_name}
        COMMAND ${CMAKE_SOURCE_DIR}/test/${test_name}.sh ${CMAKE_SOURCE_DIR}
    )
endmacro()

# Macro for MPI tests requiring runtime MPIRUN precheck.
# mode: min:<N> or exact:<N>
macro(add_hphi_mpi_test test_name mode)
    add_test(
        NAME ${test_name}
        COMMAND /bin/sh ${CMAKE_SOURCE_DIR}/test/run_with_mpi_precheck.sh
                ${mode} ${CMAKE_SOURCE_DIR}/test/${test_name}.sh
    )
    set_tests_properties(
        ${test_name}
        PROPERTIES SKIP_RETURN_CODE 77
    )
endmacro()

# Macro for MPI tests with CMAKE_SOURCE_DIR argument and precheck.
macro(add_hphi_mpi_test_with_srcdir test_name mode)
    add_test(
        NAME ${test_name}
        COMMAND /bin/sh ${CMAKE_SOURCE_DIR}/test/run_with_mpi_precheck.sh
                ${mode} ${CMAKE_SOURCE_DIR}/test/${test_name}.sh ${CMAKE_SOURCE_DIR}
    )
    set_tests_properties(
        ${test_name}
        PROPERTIES SKIP_RETURN_CODE 77
    )
endmacro()

#
# Lanczos
#
add_hphi_test(lanczos_hubbard_square)
add_hphi_test(lanczos_hubbard_square_calchs2)
add_hphi_test(lanczos_hubbard_square_ncond_calchs2)
add_hphi_test(lanczos_hubbard_chain_polarized_calchs2)
add_hphi_test(calchs2_unsupported_reject)
add_hphi_test(lanczos_kondo_chain)
add_hphi_test(lanczos_spin_kagome)
add_hphi_test(lanczos_spin_ladder_DM)
add_hphi_test(lanczos_genspin_ladder)
add_hphi_test(lanczos_hubbardgc_tri)
add_hphi_test(lanczos_kondogc_chain)
add_hphi_test(lanczos_spingc_honey)
add_hphi_test(lanczos_spingc_hcor)
add_hphi_test(lanczos_genspingc_ladder)
add_hphi_test(lanczos_tjgc_dimension)
add_hphi_test(lanczos_tjgc_exchange)
add_hphi_test_with_srcdir(lanczos_spinless)
add_hphi_test_with_srcdir(lanczos_spinless_GC)
add_hphi_test_with_srcdir(spinless_onebody_sigma_validation)
add_hphi_test_with_srcdir(spinless_onebody_offdiag_compare)
add_hphi_test(lanczos_hubbard_square_restart)
add_hphi_test(lanczos_kondo_chain_restart)
add_hphi_test(lanczos_spin_kagome_restart)
add_hphi_test(lanczos_genspin_ladder_restart)

# Lanczos test labels
set_tests_properties(
    lanczos_hubbard_square lanczos_hubbardgc_tri
    lanczos_hubbard_square_restart
    lanczos_hubbard_square_calchs2 lanczos_hubbard_square_ncond_calchs2
    lanczos_hubbard_chain_polarized_calchs2
    PROPERTIES LABELS "lanczos;hubbard")
set_tests_properties(
    calchs2_unsupported_reject
    PROPERTIES LABELS "lanczos;hubbard;validation")
set_tests_properties(
    lanczos_kondo_chain lanczos_kondogc_chain
    lanczos_kondo_chain_restart
    PROPERTIES LABELS "lanczos;kondo")
set_tests_properties(
    lanczos_spin_kagome lanczos_spin_ladder_DM
    lanczos_spingc_honey lanczos_spingc_hcor
    lanczos_spin_kagome_restart
    PROPERTIES LABELS "lanczos;spin")
set_tests_properties(
    lanczos_genspin_ladder lanczos_genspingc_ladder
    lanczos_genspin_ladder_restart
    PROPERTIES LABELS "lanczos;genspin")
set_tests_properties(
    lanczos_tjgc_dimension lanczos_tjgc_exchange
    PROPERTIES LABELS "lanczos;tjgc")
set_tests_properties(
    lanczos_spinless lanczos_spinless_GC
    PROPERTIES LABELS "lanczos;spinless")
set_tests_properties(
    spinless_onebody_sigma_validation
    PROPERTIES LABELS "spinless;validation")
set_tests_properties(
    spinless_onebody_offdiag_compare
    PROPERTIES LABELS "spinless;validation;onebody")

#
# LOBCG
#
add_hphi_test(lobcg_hubbard_square)
add_hphi_test(lobcg_hubbard_szunconserv)
add_hphi_test(lobcg_kondo_chain)
add_hphi_test(lobcg_spin_kagome)
add_hphi_test(lobcg_genspin_ladder)
add_hphi_test(lobcg_hubbardgc_tri)
add_hphi_test(lobcg_kondogc_chain)
add_hphi_test(lobcg_spingc_honey)
add_hphi_test(lobcg_spingc_Sy)
add_hphi_test(lobcg_spingc_hcor)
add_hphi_test(lobcg_spingc_6body)
add_hphi_test(lobcg_genspingc_ladder)
add_hphi_test_with_srcdir(lobcg_spinless)
add_hphi_test_with_srcdir(lobcg_spinless_GC)
add_hphi_test(lobcg_hubbard_square_restart)
add_hphi_test(lobcg_kondo_chain_restart)
add_hphi_test(lobcg_spin_kagome_restart)
add_hphi_test(lobcg_genspin_ladder_restart)

# LOBCG test labels
set_tests_properties(
    lobcg_hubbard_square lobcg_hubbard_szunconserv lobcg_hubbardgc_tri
    lobcg_hubbard_square_restart
    PROPERTIES LABELS "lobcg;hubbard")
set_tests_properties(
    lobcg_kondo_chain lobcg_kondogc_chain
    lobcg_kondo_chain_restart
    PROPERTIES LABELS "lobcg;kondo")
set_tests_properties(
    lobcg_spin_kagome lobcg_spingc_honey lobcg_spingc_Sy
    lobcg_spingc_hcor lobcg_spingc_6body
    lobcg_spin_kagome_restart
    PROPERTIES LABELS "lobcg;spin")
set_tests_properties(
    lobcg_genspin_ladder lobcg_genspingc_ladder
    lobcg_genspin_ladder_restart
    PROPERTIES LABELS "lobcg;genspin")
set_tests_properties(
    lobcg_spinless lobcg_spinless_GC
    PROPERTIES LABELS "lobcg;spinless")

#
# Full Diag
#
add_hphi_test(fulldiag_hubbard_chain)
add_test(
    NAME fulldiag_Hamiltonian_IO
    COMMAND ${CMAKE_SOURCE_DIR}/test/fulldiag_ham_io.sh
)
add_hphi_test(fulldiag_kondo_chain)
add_hphi_test(fulldiag_spin_tri)
add_hphi_test(fulldiag_genspin_ladder)
add_hphi_test(fulldiag_hubbardgc_tri)
add_hphi_test(fulldiag_tj_calchs1)
add_hphi_test(fulldiag_kondogc_chain)
add_hphi_test(fulldiag_spingc_tri)
add_hphi_test(fulldiag_genspingc_ladder)

# Full diag test labels
set_tests_properties(
    fulldiag_hubbard_chain fulldiag_hubbardgc_tri fulldiag_Hamiltonian_IO
    PROPERTIES LABELS "fulldiag;hubbard")
set_tests_properties(
    fulldiag_kondo_chain fulldiag_kondogc_chain
    PROPERTIES LABELS "fulldiag;kondo")
set_tests_properties(
    fulldiag_tj_calchs1
    PROPERTIES LABELS "fulldiag;tj")
set_tests_properties(
    fulldiag_spin_tri fulldiag_spingc_tri
    PROPERTIES LABELS "fulldiag;spin")
set_tests_properties(
    fulldiag_genspin_ladder fulldiag_genspingc_ladder
    PROPERTIES LABELS "fulldiag;genspin")

#
# Spectrum
#
add_hphi_test(spectrum_hubbard_square)
add_hphi_test(spectrum_kondo_chain)
add_hphi_test(spectrum_spin_kagome)
add_hphi_test(spectrum_genspin_ladder)
add_hphi_test(spectrum_hubbardgc_tri)
add_hphi_test(spectrum_kondogc_chain)
add_hphi_test(spectrum_spingc_honey)
add_hphi_test_with_srcdir(spectrum_spinless_chain)

# Spectrum test labels
set_tests_properties(
    spectrum_hubbard_square spectrum_hubbardgc_tri
    PROPERTIES LABELS "spectrum;hubbard")
set_tests_properties(
    spectrum_kondo_chain spectrum_kondogc_chain
    PROPERTIES LABELS "spectrum;kondo")
set_tests_properties(
    spectrum_spin_kagome spectrum_spingc_honey
    PROPERTIES LABELS "spectrum;spin")
set_tests_properties(
    spectrum_genspin_ladder
    PROPERTIES LABELS "spectrum;genspin")

#
# Time evolution
# Note: te_dc, te_pulse, te_quench depend on te_ac output directory
#
add_hphi_test(te_ac_hubbard_square)
add_hphi_test(te_dc_hubbard_square)
set_tests_properties(te_dc_hubbard_square PROPERTIES DEPENDS te_ac_hubbard_square)
add_hphi_test(te_pulse_hubbard_square)
set_tests_properties(te_pulse_hubbard_square PROPERTIES DEPENDS te_ac_hubbard_square)
add_hphi_test(te_quench_hubbard_square)
set_tests_properties(te_quench_hubbard_square PROPERTIES DEPENDS te_ac_hubbard_square)
add_hphi_test_with_srcdir(te_hubbard_chain_interall)
add_hphi_test_with_srcdir(te_hubbard_chain_interall_diagonal)
add_hphi_test_with_srcdir(te_spin_chain_interall)
add_hphi_test_with_srcdir(te_kondo_chain_interall)
add_hphi_test_with_srcdir(te_small_lanczos_max)

# Time evolution test labels
set_tests_properties(
    te_ac_hubbard_square te_dc_hubbard_square
    te_pulse_hubbard_square te_quench_hubbard_square
    te_hubbard_chain_interall te_hubbard_chain_interall_diagonal
    PROPERTIES LABELS "te;hubbard")
set_tests_properties(
    te_spin_chain_interall
    PROPERTIES LABELS "te;spin")
set_tests_properties(
    te_kondo_chain_interall
    PROPERTIES LABELS "te;kondo")
set_tests_properties(
    te_small_lanczos_max
    PROPERTIES LABELS "te;validation")

#
# mTPQ
#
add_hphi_test_with_srcdir(tpq_spin_kagome)
add_hphi_test_with_srcdir(tpq_spin_kagome_randomsphere)
add_hphi_test(tpq_small_lanczos_max)

set_tests_properties(
    tpq_spin_kagome tpq_spin_kagome_randomsphere
    PROPERTIES LABELS "tpq;spin")
set_tests_properties(
    tpq_small_lanczos_max
    PROPERTIES LABELS "tpq;validation")

#
# cTPQ
#
add_hphi_test_with_srcdir(ctpq_spin_kagome)
add_hphi_test_with_srcdir(ctpq_spin_kagome_randomsphere)
add_hphi_test_with_srcdir(read_ctpq_spin_chain)

set_tests_properties(
    ctpq_spin_kagome ctpq_spin_kagome_randomsphere read_ctpq_spin_chain
    PROPERTIES LABELS "ctpq;spin")

#
# Three body for spin 1 Heisenberg chain
#
add_hphi_test_with_srcdir(lobcg_spinone_chain_threebody)
set_tests_properties(
    lobcg_spinone_chain_threebody
    PROPERTIES LABELS "lobcg;threebody")

#
# MPI consistency tests (verify MPI and non-MPI results match)
# These tests require MPI and compare results between MPI and non-MPI runs
# Run with: MPIRUN="mpirun -np 4" ctest -R mpi_consistency
#
if(MPI_FOUND)
    # Reason markers shown in "The following tests did not run:".
    # They are intentionally skipped when MPI preconditions are not met.
    add_test(
        NAME reason_mpi_consistency_requires_MPIRUN_and_np_ge_2
        COMMAND /bin/sh ${CMAKE_SOURCE_DIR}/test/run_with_mpi_precheck.sh
                min:2 ${CMAKE_COMMAND} -E true
    )
    set_tests_properties(
        reason_mpi_consistency_requires_MPIRUN_and_np_ge_2
        PROPERTIES
        SKIP_RETURN_CODE 77
        LABELS "mpi;consistency"
    )
    add_test(
        NAME reason_mpi_batching_requires_MPIRUN_and_np_eq_16
        COMMAND /bin/sh ${CMAKE_SOURCE_DIR}/test/run_with_mpi_precheck.sh
                exact:16 ${CMAKE_COMMAND} -E true
    )
    set_tests_properties(
        reason_mpi_batching_requires_MPIRUN_and_np_eq_16
        PROPERTIES
        SKIP_RETURN_CODE 77
        LABELS "mpi;batching"
    )

    # MPI batching regression tests (MPIdouble path; script skips unless -np 16)
    add_hphi_mpi_test(lanczos_hubbard_mpidouble exact:16)
    add_hphi_mpi_test(lanczos_hubbardgc_mpidouble exact:16)

    add_hphi_mpi_test(mpi_consistency_hubbard min:2)
    add_hphi_mpi_test(mpi_consistency_hubbardgc min:2)
    add_hphi_mpi_test(mpi_consistency_spin min:2)
    add_hphi_mpi_test(mpi_consistency_spingc min:2)
    add_hphi_mpi_test(mpi_consistency_kondo min:2)
    add_hphi_mpi_test(mpi_consistency_te_hubbard min:2)
    add_hphi_mpi_test(mpi_consistency_te_spin min:2)
    add_hphi_mpi_test_with_srcdir(mpi_consistency_spinless min:2)
    add_hphi_mpi_test_with_srcdir(mpi_consistency_spinless_GC min:2)
    add_hphi_mpi_test(mpi_consistency_te_kondo min:2)
    add_hphi_mpi_test(mpi_consistency_te_hubbardgc_twobody min:4)
    add_hphi_mpi_test(mpi_consistency_te_hubbard_twobody min:4)
    add_hphi_mpi_test(mpi_nobatch_equivalence min:4)

    # MPI consistency test labels
    set_tests_properties(
        mpi_nobatch_equivalence
        PROPERTIES LABELS "mpi;nobatch")
    set_tests_properties(
        mpi_consistency_hubbard mpi_consistency_hubbardgc mpi_consistency_te_hubbard
        mpi_consistency_te_hubbardgc_twobody mpi_consistency_te_hubbard_twobody
        PROPERTIES LABELS "mpi;consistency;hubbard")
    set_tests_properties(
        mpi_consistency_spin mpi_consistency_spingc mpi_consistency_te_spin
        PROPERTIES LABELS "mpi;consistency;spin")
    set_tests_properties(
        mpi_consistency_kondo mpi_consistency_te_kondo
        PROPERTIES LABELS "mpi;consistency;kondo")
    set_tests_properties(
        mpi_consistency_spinless mpi_consistency_spinless_GC
        PROPERTIES LABELS "mpi;consistency;spinless")
    set_tests_properties(
        lanczos_hubbard_mpidouble lanczos_hubbardgc_mpidouble
        PROPERTIES LABELS "mpi;lanczos;hubbard;batching")
endif(MPI_FOUND)
