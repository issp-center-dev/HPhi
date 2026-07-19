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

add_executable(unittest_symmetry_basis
    ${CMAKE_SOURCE_DIR}/test/unittest_symmetry_basis.c
    ${CMAKE_SOURCE_DIR}/src/symmetry_basis.c
    ${CMAKE_SOURCE_DIR}/src/symmetry_matvec_plan.c
    ${CMAKE_SOURCE_DIR}/src/bitcalc.c
    ${CMAKE_SOURCE_DIR}/src/ErrorMessage.c)
target_include_directories(unittest_symmetry_basis
    PRIVATE ${CMAKE_SOURCE_DIR}/src/include)
target_compile_definitions(unittest_symmetry_basis PRIVATE
    DSFMT_MEXP=19937 HPHI_SYMMETRY_PARITY_VERIFY)
target_link_libraries(unittest_symmetry_basis m)
if(MPI_FOUND)
    target_link_libraries(unittest_symmetry_basis ${MPI_C_LIBRARIES})
endif()
if(TARGET OpenMP::OpenMP_C)
    target_link_libraries(unittest_symmetry_basis OpenMP::OpenMP_C)
endif()
add_test(NAME unittest_symmetry_basis
    COMMAND unittest_symmetry_basis)
set_tests_properties(unittest_symmetry_basis
    PROPERTIES LABELS "unit;symmetry")

#
# Lanczos
#
add_hphi_test(lanczos_hubbard_square)
add_hphi_test(lanczos_hubbard_square_calchs2)
add_hphi_test(lanczos_hubbard_square_ncond_calchs2)
add_hphi_test(modpara_parse_validation)
set_tests_properties(modpara_parse_validation
    PROPERTIES LABELS "validation;input")
add_hphi_test(lanczos_hubbard_chain_polarized_calchs2)
add_hphi_test(calchs2_unsupported_reject)
add_hphi_test(threebody_calcmodel_validation)
add_hphi_test(lanczos_small_dimension_regression)
add_hphi_test(symmetry_spin_chain_validation)
add_hphi_test(symmetry_spin_chain_standard_generation)
add_hphi_test(symmetry_spinless_fermion_lanczos)
add_hphi_test(symmetry_spinless_fermion_standard_generation)
add_hphi_test(symmetry_hubbard_lanczos)
add_hphi_test(symmetry_hubbard_standard_generation)
add_hphi_test(lanczos_kondo_chain)
add_hphi_test(lanczos_spin_kagome)
add_hphi_test(lanczos_spin_calchs_minus1)
add_hphi_test(symmetry_spin_chain_lanczos)
add_hphi_test(lanczos_spin_ladder_DM)
add_hphi_test(lanczos_genspin_ladder)
add_hphi_test(lanczos_hubbardgc_tri)
add_hphi_test(lanczos_hubbardgc_nbody_interall)
add_hphi_test(lanczos_hubbardgc_nbodyg)
add_hphi_test(lanczos_hubbardgc_anomalous_pair)
add_hphi_test(lanczos_hubbard_nbody_interall)
add_hphi_test(lanczos_hubbard_nbodyg)
add_hphi_test(lanczos_kondogc_chain)
add_hphi_test(lanczos_tj_kondo_nbody_interall)
add_hphi_test(lanczos_tj_nconserved_nbodyg)
add_hphi_test(lanczos_spingc_honey)
add_hphi_test(lanczos_spingc_hcor)
add_hphi_test(lanczos_spingc_nbody_interall)
add_hphi_test(lanczos_spingc_nbody_interall_n1)
add_hphi_test(lanczos_spingc_nbodyg)
add_hphi_test(lanczos_spin_nbody_interall)
add_hphi_test(lanczos_spingc_spinone_nbody_interall)
add_hphi_test(lanczos_spingc_spinone_nbodyg)
add_hphi_test(lanczos_spin_spinone_nbody_interall)
add_hphi_test(lanczos_spin_spinone_nbodyg)
add_hphi_test(fulldiag_spin_spinone_nbodyg)
add_hphi_test(nbody_spingc_spinone_validation)
add_hphi_test(lanczos_genspingc_ladder)
add_hphi_test(lanczos_tjgc_dimension)
add_hphi_test(lanczos_tjgc_exchange)
add_hphi_test(nbody_interall_validation)
add_hphi_test(nbodyg_validation)
add_hphi_test(nbody_hubbardgc_validation)
add_hphi_test(anomalous_hubbardgc_validation)
add_hphi_test(nbody_hubbard_validation)
add_hphi_test(nbody_tj_kondo_validation)
add_hphi_test_with_srcdir(lanczos_spinless)
add_hphi_test_with_srcdir(lanczos_spinless_GC)
add_hphi_test_with_srcdir(lanczos_spinless_nbody_interall)
add_hphi_test_with_srcdir(lanczos_spinless_nbodyg)
add_hphi_test_with_srcdir(nbody_spinless_validation)
add_hphi_test_with_srcdir(spinless_interall_diagonal)
add_hphi_test_with_srcdir(lanczos_spinless_calchs0)
add_hphi_test_with_srcdir(spinless_onebody_sigma_validation)
add_hphi_test_with_srcdir(spinless_onebody_offdiag_compare)
add_hphi_test_with_srcdir(green_output_format)
add_hphi_test(lanczos_hubbard_square_restart)
add_hphi_test(lanczos_kondo_chain_restart)
add_hphi_test(lanczos_spin_kagome_restart)
add_hphi_test(lanczos_genspin_ladder_restart)

# Lanczos test labels
set_tests_properties(
    lanczos_hubbard_square lanczos_hubbardgc_tri
    lanczos_hubbardgc_nbody_interall lanczos_hubbardgc_nbodyg
    lanczos_hubbard_nbody_interall lanczos_hubbard_nbodyg
    lanczos_hubbard_square_restart
    lanczos_hubbard_square_calchs2 lanczos_hubbard_square_ncond_calchs2
    lanczos_hubbard_chain_polarized_calchs2
    PROPERTIES LABELS "lanczos;hubbard")
set_tests_properties(
    lanczos_hubbardgc_anomalous_pair
    PROPERTIES LABELS "lanczos;fulldiag;te;hubbard;anomalous")
set_tests_properties(
    calchs2_unsupported_reject
    PROPERTIES LABELS "lanczos;hubbard;validation")
set_tests_properties(
    threebody_calcmodel_validation
    PROPERTIES LABELS "lanczos;spin;validation")
set_tests_properties(
    lanczos_small_dimension_regression
    PROPERTIES LABELS "lanczos;spin;validation")
set_tests_properties(
    symmetry_spin_chain_validation
    PROPERTIES LABELS "symmetry;spin;validation")
set_tests_properties(
    symmetry_spin_chain_standard_generation
    PROPERTIES LABELS "symmetry;spin;standard")
set_tests_properties(
    symmetry_spinless_fermion_lanczos
    PROPERTIES LABELS "symmetry;spinless;lanczos")
set_tests_properties(
    symmetry_spinless_fermion_standard_generation
    PROPERTIES LABELS "symmetry;spinless;standard")
set_tests_properties(
    symmetry_hubbard_lanczos
    PROPERTIES LABELS "symmetry;hubbard;lanczos;validation")
set_tests_properties(
    symmetry_hubbard_standard_generation
    PROPERTIES LABELS "symmetry;hubbard;standard")
set_tests_properties(
    lanczos_kondo_chain lanczos_kondogc_chain
    lanczos_kondo_chain_restart
    PROPERTIES LABELS "lanczos;kondo")
set_tests_properties(
    lanczos_spin_kagome lanczos_spin_ladder_DM
    lanczos_spin_calchs_minus1
    lanczos_spingc_honey lanczos_spingc_hcor
    lanczos_spingc_nbodyg
    lanczos_spin_nbody_interall
    lanczos_spin_kagome_restart
    PROPERTIES LABELS "lanczos;spin")
set_tests_properties(
    symmetry_spin_chain_lanczos
    PROPERTIES LABELS "symmetry;spin;lanczos")
set_tests_properties(
    lanczos_genspin_ladder lanczos_genspingc_ladder
    lanczos_genspin_ladder_restart
    PROPERTIES LABELS "lanczos;genspin")
set_tests_properties(
    lanczos_tjgc_dimension lanczos_tjgc_exchange
    PROPERTIES LABELS "lanczos;tjgc")
set_tests_properties(
    lanczos_tj_kondo_nbody_interall
    PROPERTIES LABELS "lanczos;tj;tjgc;kondo;nbody")
set_tests_properties(
    lanczos_tj_nconserved_nbodyg
    PROPERTIES LABELS "lanczos;tj;nbody")
set_tests_properties(
    lanczos_spinless lanczos_spinless_GC lanczos_spinless_calchs0
    PROPERTIES LABELS "lanczos;spinless")
set_tests_properties(
    lanczos_spinless_nbody_interall lanczos_spinless_nbodyg
    PROPERTIES LABELS "lanczos;spinless;nbody")
set_tests_properties(
    nbodyg_validation
    PROPERTIES LABELS "spin;validation;nbody")
set_tests_properties(
    nbody_spingc_spinone_validation
    PROPERTIES LABELS "spinone;genspin;validation;nbody")
set_tests_properties(
    nbody_hubbardgc_validation nbody_hubbard_validation
    PROPERTIES LABELS "hubbard;validation;nbody")
set_tests_properties(
    anomalous_hubbardgc_validation
    PROPERTIES LABELS "hubbard;validation;anomalous")
set_tests_properties(
    nbody_tj_kondo_validation
    PROPERTIES LABELS "tj;tjgc;kondo;validation;nbody")
set_tests_properties(
    nbody_spinless_validation
    PROPERTIES LABELS "spinless;validation;nbody")
set_tests_properties(
    lanczos_spingc_spinone_nbody_interall
    lanczos_spin_spinone_nbody_interall
    PROPERTIES LABELS "lanczos;genspin;spinone;nbody")
set_tests_properties(
    lanczos_spingc_spinone_nbodyg
    lanczos_spin_spinone_nbodyg
    PROPERTIES LABELS "lanczos;genspin;spinone;nbody")
set_tests_properties(
    fulldiag_spin_spinone_nbodyg
    PROPERTIES LABELS "fulldiag;genspin;spinone;nbody")
set_tests_properties(
    spinless_onebody_sigma_validation
    PROPERTIES LABELS "spinless;validation")
set_tests_properties(
    spinless_onebody_offdiag_compare
    PROPERTIES LABELS "spinless;validation;onebody")
set_tests_properties(
    green_output_format
    PROPERTIES LABELS "tpq;te;fulldiag;green;validation")

#
# LOBCG
#
add_hphi_test(lobcg_hubbard_square)
add_hphi_test(lobcg_hubbard_szunconserv)
add_hphi_test(lobcg_kondo_chain)
add_hphi_test(lobcg_spin_kagome)
add_hphi_test(symmetry_spin_chain_lobcg)
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
    symmetry_spin_chain_lobcg
    PROPERTIES LABELS "symmetry;spin;lobcg")
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
add_hphi_test(kondo_nconserved_dimension)
add_hphi_test(fulldiag_spin_tri)
add_hphi_test(fulldiag_genspin_ladder)
add_hphi_test(fulldiag_hubbardgc_tri)
add_hphi_test(fulldiag_hubbardgc_nbody_interall)
add_hphi_test(fulldiag_hubbard_nbody_interall)
add_hphi_test(fulldiag_tj_calchs1)
add_hphi_test(fulldiag_tj_kondo_nbody_interall)
add_hphi_test(fulldiag_kondogc_chain)
add_hphi_test(fulldiag_spingc_tri)
add_hphi_test(fulldiag_spingc_nbody_interall_compare)
add_hphi_test(fulldiag_genspingc_ladder)

# Full diag test labels
set_tests_properties(
    fulldiag_hubbard_chain fulldiag_hubbardgc_tri fulldiag_Hamiltonian_IO
    fulldiag_hubbardgc_nbody_interall
    fulldiag_hubbard_nbody_interall
    PROPERTIES LABELS "fulldiag;hubbard")
set_tests_properties(
    fulldiag_kondo_chain fulldiag_kondogc_chain
    PROPERTIES LABELS "fulldiag;kondo")
set_tests_properties(
    kondo_nconserved_dimension
    PROPERTIES LABELS "fulldiag;kondo;validation")
set_tests_properties(
    fulldiag_tj_calchs1
    PROPERTIES LABELS "fulldiag;tj")
set_tests_properties(
    fulldiag_tj_kondo_nbody_interall
    PROPERTIES LABELS "fulldiag;tj;tjgc;kondo;nbody")
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

# Off-diagonal dynamical Green function (bra excitation)
add_hphi_test(spectrum_hubbard_offdiag_single)
add_hphi_test(spectrum_hubbard_offdiag_pair_density)
add_hphi_test(spectrum_spin_offdiag_szsz)
add_hphi_test(spectrum_spingc_offdiag_szsz)
add_hphi_test(spectrum_offdiag_invalid_reject)
add_hphi_test(spectrum_bicg_nonconvergence_reject)
add_hphi_test(spectrum_loop_multiop_reader_validation)
add_hphi_test(spectrum_hubbard_nconserved_calchs2_boundary)

# Spectrum test labels
set_tests_properties(
    spectrum_hubbard_square spectrum_hubbardgc_tri
    PROPERTIES LABELS "spectrum;hubbard")
set_tests_properties(
    spectrum_hubbard_offdiag_single spectrum_hubbard_offdiag_pair_density
    PROPERTIES LABELS "spectrum;hubbard;offdiag")
set_tests_properties(
    spectrum_spin_offdiag_szsz spectrum_spingc_offdiag_szsz
    PROPERTIES LABELS "spectrum;spin;offdiag")
set_tests_properties(
    spectrum_offdiag_invalid_reject
    PROPERTIES LABELS "spectrum;offdiag;validation")
set_tests_properties(
    spectrum_bicg_nonconvergence_reject
    PROPERTIES LABELS "spectrum;hubbard;offdiag;validation")
set_tests_properties(
    spectrum_loop_multiop_reader_validation
    PROPERTIES LABELS "spectrum;hubbard;validation")
set_tests_properties(
    spectrum_hubbard_nconserved_calchs2_boundary
    PROPERTIES LABELS "spectrum;hubbard;offdiag;validation")

# Unit test: off-diagonal excitation sector-shift truth table
add_executable(unittest_excitation_sector_shift
    ${CMAKE_SOURCE_DIR}/test/unittest_excitation_sector_shift.c
    ${CMAKE_SOURCE_DIR}/src/ExcitationSectorShift.c)
target_include_directories(unittest_excitation_sector_shift
    PRIVATE ${CMAKE_SOURCE_DIR}/src/include)
target_compile_definitions(unittest_excitation_sector_shift PRIVATE DSFMT_MEXP=19937)
add_test(NAME unittest_excitation_sector_shift
    COMMAND unittest_excitation_sector_shift)
set_tests_properties(unittest_excitation_sector_shift
    PROPERTIES LABELS "unit;offdiag")

# Unit test: zero-size 2d/3d setmemory allocations
add_executable(unittest_setmemory_zero_size_common
    ${CMAKE_SOURCE_DIR}/test/unittest_setmemory_zero_size.c
    ${CMAKE_SOURCE_DIR}/src/common/setmemory.c)
target_include_directories(unittest_setmemory_zero_size_common
    PRIVATE ${CMAKE_SOURCE_DIR}/src/common)
add_test(NAME unittest_setmemory_zero_size_common
    COMMAND unittest_setmemory_zero_size_common)
set_tests_properties(unittest_setmemory_zero_size_common
    PROPERTIES LABELS "unit;setmemory")

# Unit test: komega BiCG non-finite guard
add_executable(unittest_komega_bicg_nonfinite
    ${CMAKE_SOURCE_DIR}/test/unittest_komega_bicg_nonfinite.c)
target_include_directories(unittest_komega_bicg_nonfinite
    PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(unittest_komega_bicg_nonfinite komega)
if(MPI_FOUND)
    target_link_libraries(unittest_komega_bicg_nonfinite ${MPI_C_LIBRARIES})
endif(MPI_FOUND)
add_test(NAME unittest_komega_bicg_nonfinite
    COMMAND unittest_komega_bicg_nonfinite)
set_tests_properties(unittest_komega_bicg_nonfinite
    PROPERTIES LABELS "unit;spectrum")

# MPI consistency for the off-diagonal spectrum (4-site Hubbard needs np=4)
add_hphi_mpi_test(spectrum_hubbard_offdiag_single_mpi exact:4)
set_tests_properties(spectrum_hubbard_offdiag_single_mpi
    PROPERTIES LABELS "spectrum;hubbard;offdiag;mpi")
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
add_hphi_test_with_srcdir(te_spinless_diagonal_reject)

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
set_tests_properties(
    te_spinless_diagonal_reject
    PROPERTIES LABELS "te;spinless;validation")

#
# mTPQ
#
add_hphi_test_with_srcdir(tpq_spin_kagome)
add_hphi_test_with_srcdir(tpq_spin_kagome_randomsphere)
add_hphi_test(tpq_small_lanczos_max)
add_hphi_test(tpq_restart_missing_vector_fallback)

set_tests_properties(
    tpq_spin_kagome tpq_spin_kagome_randomsphere
    PROPERTIES LABELS "tpq;spin")
set_tests_properties(
    tpq_small_lanczos_max
    PROPERTIES LABELS "tpq;validation")
set_tests_properties(
    tpq_restart_missing_vector_fallback
    PROPERTIES LABELS "tpq;ctpq;validation")

#
# cTPQ
#
add_hphi_test_with_srcdir(ctpq_spin_kagome)
add_hphi_test_with_srcdir(ctpq_spin_kagome_randomsphere)
add_hphi_test_with_srcdir(read_ctpq_spin_chain)
add_hphi_test(ctpq_invtemp_validation)

set_tests_properties(
    ctpq_spin_kagome ctpq_spin_kagome_randomsphere read_ctpq_spin_chain
    PROPERTIES LABELS "ctpq;spin")
set_tests_properties(
    ctpq_invtemp_validation
    PROPERTIES LABELS "ctpq;validation")

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
    add_hphi_mpi_test(mpi_consistency_kondo_nconserved_cg exact:4)
    add_hphi_mpi_test(mpi_consistency_te_hubbard min:2)
    add_hphi_mpi_test(mpi_consistency_te_spin min:2)
    add_hphi_mpi_test_with_srcdir(mpi_consistency_spinless min:2)
    add_hphi_mpi_test_with_srcdir(mpi_consistency_spinless_GC min:2)
    add_hphi_mpi_test(mpi_consistency_te_kondo min:2)
    add_hphi_mpi_test_with_srcdir(mpi_consistency_threebody_generalspin exact:3)
    add_hphi_mpi_test(mpi_consistency_te_hubbardgc_twobody min:4)
    add_hphi_mpi_test(mpi_consistency_te_hubbard_twobody min:4)
    add_hphi_mpi_test(mpi_nobatch_equivalence min:4)
    # batched-vs-nobatch equivalence for the paths the base test does not reach:
    # MPIdouble (np=16), HubbardGC off-diagonal InterAll, SpinGC PairLift, Spinless.
    add_hphi_mpi_test(mpi_nobatch_equivalence_mpidouble exact:16)
    add_hphi_mpi_test(mpi_nobatch_equivalence_interall exact:4)
    add_hphi_mpi_test(mpi_nobatch_equivalence_interall_selfterm exact:4)
    add_hphi_mpi_test(mpi_nobatch_equivalence_pairlift exact:4)
    add_hphi_mpi_test_with_srcdir(mpi_nobatch_equivalence_spinless exact:4)
    add_hphi_mpi_test(mpi_nbody_interall_spingc exact:4)
    add_hphi_mpi_test(mpi_nbodyg_spingc exact:4)
    add_hphi_mpi_test(mpi_nbody_interall_hubbardgc exact:4)
    add_hphi_mpi_test(mpi_nbodyg_hubbardgc exact:4)
    add_hphi_mpi_test(mpi_anomalous_hubbardgc exact:4)
    add_hphi_mpi_test(mpi_nbody_interall_hubbard exact:4)
    add_hphi_mpi_test(mpi_nbodyg_hubbard exact:4)
    add_hphi_mpi_test(mpi_nbody_interall_kondo_nconserved exact:4)
    add_hphi_mpi_test_with_srcdir(mpi_nbody_interall_spinless exact:4)
    add_hphi_mpi_test_with_srcdir(mpi_nbodyg_spinless exact:4)
    add_hphi_mpi_test(mpi_nbody_interall_spin exact:4)
    add_hphi_mpi_test(mpi_nbodyg_spin exact:4)
    add_hphi_mpi_test(mpi_nbody_interall_spingc_spinone exact:3)
    add_hphi_mpi_test(mpi_nbody_interall_spingc_spinone_np9 exact:9)
    add_hphi_mpi_test(mpi_nbodyg_spingc_spinone exact:3)
    add_hphi_mpi_test(mpi_nbodyg_spingc_spinone_np9 exact:9)
    add_hphi_mpi_test(mpi_nbody_interall_spin_spinone exact:3)
    add_hphi_mpi_test(mpi_nbody_interall_spin_spinone_np9 exact:9)
    add_hphi_mpi_test(mpi_nbodyg_spin_spinone exact:3)
    add_hphi_mpi_test(mpi_nbodyg_spin_spinone_np9 exact:9)

    # MPI consistency test labels
    set_tests_properties(
        mpi_nobatch_equivalence
        mpi_nobatch_equivalence_mpidouble
        mpi_nobatch_equivalence_interall
        mpi_nobatch_equivalence_interall_selfterm
        mpi_nobatch_equivalence_pairlift
        mpi_nobatch_equivalence_spinless
        PROPERTIES LABELS "mpi;nobatch")
    set_tests_properties(
        mpi_consistency_hubbard mpi_consistency_hubbardgc mpi_consistency_te_hubbard
        mpi_consistency_te_hubbardgc_twobody mpi_consistency_te_hubbard_twobody
        mpi_nbody_interall_hubbardgc mpi_nbodyg_hubbardgc
        mpi_anomalous_hubbardgc
        mpi_nbody_interall_hubbard mpi_nbodyg_hubbard
        PROPERTIES LABELS "mpi;consistency;hubbard")
    set_tests_properties(
        mpi_consistency_spin mpi_consistency_spingc mpi_consistency_te_spin
        mpi_consistency_threebody_generalspin
        PROPERTIES LABELS "mpi;consistency;spin")
    set_tests_properties(
        mpi_nbody_interall_spingc mpi_nbodyg_spingc
        mpi_nbody_interall_spin mpi_nbodyg_spin
        PROPERTIES LABELS "mpi;consistency;spin;nbody")
    set_tests_properties(
        mpi_nbody_interall_spingc_spinone
        mpi_nbody_interall_spingc_spinone_np9
        mpi_nbodyg_spingc_spinone
        mpi_nbodyg_spingc_spinone_np9
        mpi_nbody_interall_spin_spinone
        mpi_nbody_interall_spin_spinone_np9
        mpi_nbodyg_spin_spinone
        mpi_nbodyg_spin_spinone_np9
        PROPERTIES LABELS "mpi;consistency;genspin;spinone;nbody")
    set_tests_properties(
        mpi_consistency_kondo mpi_consistency_kondo_nconserved_cg
        mpi_consistency_te_kondo
        PROPERTIES LABELS "mpi;consistency;kondo")
    set_tests_properties(
        mpi_nbody_interall_kondo_nconserved
        PROPERTIES LABELS "mpi;consistency;kondo;nbody")
    set_tests_properties(
        mpi_consistency_spinless mpi_consistency_spinless_GC
        PROPERTIES LABELS "mpi;consistency;spinless")
    set_tests_properties(
        mpi_nbody_interall_spinless mpi_nbodyg_spinless
        PROPERTIES LABELS "mpi;consistency;spinless;nbody")
    set_tests_properties(
        lanczos_hubbard_mpidouble lanczos_hubbardgc_mpidouble
        PROPERTIES LABELS "mpi;lanczos;hubbard;batching")
endif(MPI_FOUND)
