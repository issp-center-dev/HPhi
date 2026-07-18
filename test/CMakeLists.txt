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
add_hphi_test(threebody_calcmodel_validation)
add_hphi_test(lanczos_kondo_chain)
add_hphi_test(lanczos_spin_kagome)
add_hphi_test(lanczos_spin_calchs_minus1)
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
add_hphi_test_with_srcdir(check_expec_local_calls)
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
add_hphi_test(fulldiag_solver_keyword)
if(USE_SCALAPACK)
  set_property(TEST fulldiag_solver_keyword APPEND PROPERTY ENVIRONMENT "HPHI_HAS_SCALAPACK=1")
endif(USE_SCALAPACK)
if(USE_ELPA)
  set_property(TEST fulldiag_solver_keyword APPEND PROPERTY ENVIRONMENT "HPHI_HAS_ELPA=1")
endif(USE_ELPA)
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
add_hphi_test(fulldiag_spingc_gamma)

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
    fulldiag_spin_tri fulldiag_spingc_tri fulldiag_spingc_gamma
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

# Unit test: ExpecMode 2 trace-kernel mapping-probe extraction (phase 3b Task 2)
# and the ONEBODY streaming kernel (phase 3b Task 3). Serial, no MPI: GC
# models only (HubbardGC n=256, SpinGC-half n=16), so the bare-bit basis needs
# no list_1/list_2 and this runs as a plain executable in BOTH the MPI and
# non-MPI builds (registered here, outside if(MPI_FOUND)).
# Links only the real sources in the extraction closure -- expec_trace.c (the
# drivers, plus, since Task 3, TraceStreamOneBody() and the ONEBODY output
# phase), the two Core files (the co-located *_TraceProbe adapters + element
# functions), rearray_interactions.c (the REAL two-body reordering, its own TU
# precisely so this test does not need a drift-prone copy), bitcalc.c
# (GetSplitBitByModel/GetOffComp/SgnBit), global.c (which provides the
# list_1/Tpow globals plus `myrank` and stdoutMPI), and ErrorMessage.c (the
# cErr* strings bitcalc references). Task 3 added wrapperMPI.c/FileIO.c/
# green_output.c/splash.c to this list: expec_trace.c's
# (never-called-by-this-test) output phase now references
# childfopenMPI()/GreenOutput*() symbols, so the TU needs them at LINK time
# even though the test only exercises TraceStreamOneBody() (pure buffer
# fill, no I/O) and TraceBuildPlan(), never expec_trace_owned_states()
# itself; splash.c is transitively required because wrapperMPI.c's
# (likewise never-called) InitializeMPI() references splash(). The
# MPI-heavy element libraries are still NOT pulled in.
add_executable(expec_trace_map_check unit/expec_trace_map_check.c
    ${CMAKE_SOURCE_DIR}/src/expec_trace.c
    ${CMAKE_SOURCE_DIR}/src/mltplyHubbardCore.c
    ${CMAKE_SOURCE_DIR}/src/mltplySpinCore.c
    ${CMAKE_SOURCE_DIR}/src/rearray_interactions.c
    ${CMAKE_SOURCE_DIR}/src/bitcalc.c
    ${CMAKE_SOURCE_DIR}/src/global.c
    ${CMAKE_SOURCE_DIR}/src/ErrorMessage.c
    ${CMAKE_SOURCE_DIR}/src/wrapperMPI.c
    ${CMAKE_SOURCE_DIR}/src/FileIO.c
    ${CMAKE_SOURCE_DIR}/src/green_output.c
    ${CMAKE_SOURCE_DIR}/src/splash.c)
target_include_directories(expec_trace_map_check PRIVATE
    ${CMAKE_SOURCE_DIR}/src/include)
target_compile_definitions(expec_trace_map_check PRIVATE DSFMT_MEXP=19937)
target_link_libraries(expec_trace_map_check m)
# expec_trace.c's clock_gettime may live in librt on older glibc -- same
# check-and-link as the HPhi target (src/CMakeLists.txt).
include(CheckLibraryExists)
check_library_exists(rt clock_gettime "" HPHI_CLOCK_GETTIME_IN_RT)
if(HPHI_CLOCK_GETTIME_IN_RT)
  target_link_libraries(expec_trace_map_check rt)
endif()
# On an MPI build the top-level add_definitions(-DMPI) makes the linked-in
# src/green_output.c contain real MPI calls (GreenOutputMergePartials()'s
# MPI_Gather/MPI_Bcast), so the MPI libraries must be on this target's link
# line -- same pattern/condition as the green_partial_merge_check target
# below. The test itself stays serial and launcher-free: it never calls
# GreenOutputMergePartials() (nor expec_trace_owned_states(), whose output
# phase is the only expec_trace.c path into green_output.c), so no MPI
# function is ever reached at runtime and MPI_Init is not required.
if(MPI_FOUND)
    target_link_libraries(expec_trace_map_check ${MPI_C_LIBRARIES})
endif(MPI_FOUND)
add_test(NAME expec_trace_map_check COMMAND expec_trace_map_check)
set_tests_properties(expec_trace_map_check
    PROPERTIES LABELS "unit;fulldiag;expecmode")

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

    # green_output partial-output/merge unit test (phase 3a Task 7). Unlike
    # elpa_eigen_check/elpa_redist_check/elpa_statepanel_check above, the
    # code this test exercises (GreenOutputMergePartials() and the partial
    # session API in src/green_output.c) is guarded only by "#ifdef MPI",
    # NOT by "_SCALAPACK"/"_ELPA" -- it has nothing to do with the
    # ScaLAPACK/ELPA eigensolver path, only with MPI-rank-local partial
    # Green-output files. So it is registered here, under plain
    # if(MPI_FOUND), and actually RUNS on any MPI build (ELPA/ScaLAPACK
    # present or not) rather than being gated behind if(USE_ELPA).
    # It links the real src/green_output.c + src/FileIO.c + src/wrapperMPI.c
    # + src/global.c + src/ErrorMessage.c (no HPhi-global stub objects
    # substituted) so the test exercises the actual merge/publish code
    # path, including the real childfopenMPI()/fopenMPI() ExpecLocal
    # per-rank file I/O the Mode 1 driver relies on.
    add_executable(green_partial_merge_check unit/green_partial_merge_check.c
        ${CMAKE_SOURCE_DIR}/src/green_output.c
        ${CMAKE_SOURCE_DIR}/src/FileIO.c
        ${CMAKE_SOURCE_DIR}/src/wrapperMPI.c
        ${CMAKE_SOURCE_DIR}/src/global.c
        ${CMAKE_SOURCE_DIR}/src/ErrorMessage.c)
    target_include_directories(green_partial_merge_check PRIVATE
        ${CMAKE_SOURCE_DIR}/src/include)
    target_compile_definitions(green_partial_merge_check PRIVATE -DMPI -DGREEN_OUTPUT_TESTING)
    target_link_libraries(green_partial_merge_check ${MPI_C_LIBRARIES} m)
    # Same run_with_mpi_precheck.sh + SKIP_RETURN_CODE 77 convention as the
    # other MPI unit tests. min:2 (not exact:2): every scenario in
    # unit/green_partial_merge_check.c is written to generalize to any
    # nprocs >= 2 (rank 0 / highest rank play fixed roles).
    add_test(
        NAME green_partial_merge_check
        COMMAND /bin/sh ${CMAKE_SOURCE_DIR}/test/run_with_mpi_precheck.sh
                min:2 ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                ${MPIEXEC_PREFLAGS} $<TARGET_FILE:green_partial_merge_check>
                ${MPIEXEC_POSTFLAGS}
    )
    set_tests_properties(
        green_partial_merge_check
        PROPERTIES SKIP_RETURN_CODE 77 LABELS "fulldiag;green;unit;mpi"
    )
endif(MPI_FOUND)

if(USE_ELPA)
    configure_file(
        ${CMAKE_SOURCE_DIR}/test/cmake/elpa_configure_cache.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/elpa_configure_cache.cmake
        @ONLY)
    add_test(
        NAME elpa_configure_check
        COMMAND ${CMAKE_COMMAND}
                -DHPHI_SOURCE_DIR=${CMAKE_SOURCE_DIR}
                -DHPHI_TEST_ROOT=${CMAKE_CURRENT_BINARY_DIR}/elpa_configure_check
                -DHPHI_PRELOAD_CACHE=${CMAKE_CURRENT_BINARY_DIR}/elpa_configure_cache.cmake
                -P ${CMAKE_SOURCE_DIR}/test/cmake/elpa_configure_check.cmake)
    set_tests_properties(
        elpa_configure_check
        PROPERTIES LABELS "elpa;cmake;validation")

    add_hphi_mpi_test(fulldiag_elpa_hubbard_chain min:2)
    add_hphi_mpi_test(fulldiag_output_open_failure exact:2)
    set_tests_properties(
        fulldiag_output_open_failure
        PROPERTIES LABELS "fulldiag;elpa;validation" TIMEOUT 30)

    add_executable(elpa_eigen_check unit/elpa_eigen_check.c
        ${CMAKE_SOURCE_DIR}/src/matrixscalapack.c
        ${CMAKE_SOURCE_DIR}/src/matrixlapack_elpa.c)
    target_include_directories(elpa_eigen_check PRIVATE
        ${CMAKE_SOURCE_DIR}/src/include ${ELPA_INCLUDE_DIRS})
    target_compile_options(elpa_eigen_check PRIVATE ${ELPA_COMPILE_OPTIONS})
    target_compile_definitions(elpa_eigen_check PRIVATE -D_ELPA -D_SCALAPACK -DMPI)
    # SCALAPACK_LIBRARIES may be a single space-separated flag string from
    # config/*.cmake; normalize it into a list the same way src/CMakeLists.txt
    # does for the HPhi target, or the test target gets one malformed link arg.
    string(REGEX REPLACE "-L[ ]+" "-L" _ec_sc_libs "${SCALAPACK_LIBRARIES}")
    string(REGEX REPLACE "[ ]+" ";" _ec_sc_libs "${_ec_sc_libs}")
    target_link_libraries(elpa_eigen_check ${ELPA_LIBRARIES}
        ${_ec_sc_libs} ${LAPACK_LIBRARIES} ${MPI_C_LIBRARIES}
        ${MPI_Fortran_LIBRARIES} ${HPHI_FORTRAN_RUNTIME_LIBRARIES} m)
    # Wrapped with the same run_with_mpi_precheck.sh + SKIP_RETURN_CODE 77
    # convention as the other MPI tests in this file (see add_hphi_mpi_test),
    # rather than a bare add_test(mpiexec ...): on hosts without a usable
    # MPIRUN (or with too few ranks configured), this now skips cleanly
    # (exit 77) instead of failing/hanging trying to launch mpiexec directly.
    # Requires exactly 3 ranks (the grid assumed by unit/elpa_eigen_check.c),
    # hence mode "exact:3" like the other exact-rank-count MPI tests above.
    add_test(
        NAME elpa_eigen_check
        COMMAND /bin/sh ${CMAKE_SOURCE_DIR}/test/run_with_mpi_precheck.sh
                exact:3 ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3
                ${MPIEXEC_PREFLAGS} $<TARGET_FILE:elpa_eigen_check>
                ${MPIEXEC_POSTFLAGS}
    )
    set_tests_properties(
        elpa_eigen_check
        PROPERTIES SKIP_RETURN_CODE 77 LABELS "fulldiag;elpa;unit"
    )

    add_executable(elpa_redist_check unit/elpa_redist_check.c
        ${CMAKE_SOURCE_DIR}/src/matrixscalapack.c
        ${CMAKE_SOURCE_DIR}/src/matrixlapack_elpa.c)
    target_include_directories(elpa_redist_check PRIVATE
        ${CMAKE_SOURCE_DIR}/src/include ${ELPA_INCLUDE_DIRS})
    target_compile_options(elpa_redist_check PRIVATE ${ELPA_COMPILE_OPTIONS})
    target_compile_definitions(elpa_redist_check PRIVATE -D_ELPA -D_SCALAPACK -DMPI)
    # _ec_sc_libs is the normalized SCALAPACK_LIBRARIES list computed above
    # for elpa_eigen_check; reuse it here for the same reason.
    target_link_libraries(elpa_redist_check ${ELPA_LIBRARIES}
        ${_ec_sc_libs} ${LAPACK_LIBRARIES} ${MPI_C_LIBRARIES}
        ${MPI_Fortran_LIBRARIES} ${HPHI_FORTRAN_RUNTIME_LIBRARIES} m)
    # Wrapped with the same run_with_mpi_precheck.sh + SKIP_RETURN_CODE 77
    # convention as elpa_eigen_check above. Unlike elpa_eigen_check,
    # unit/elpa_redist_check.c builds its BLACS grid from MPI_Dims_create
    # and works with any process count, so "min:2" (rather than
    # "exact:3") is enough; launch with 2 ranks to exercise the
    # non-divisible-N redistribution path with more than one owner.
    add_test(
        NAME elpa_redist_check
        COMMAND /bin/sh ${CMAKE_SOURCE_DIR}/test/run_with_mpi_precheck.sh
                min:2 ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                ${MPIEXEC_PREFLAGS} $<TARGET_FILE:elpa_redist_check>
                ${MPIEXEC_POSTFLAGS}
    )
    set_tests_properties(
        elpa_redist_check
        PROPERTIES SKIP_RETURN_CODE 77 LABELS "fulldiag;elpa;unit"
    )

    add_executable(elpa_statepanel_check unit/elpa_statepanel_check.c
        ${CMAKE_SOURCE_DIR}/src/matrixscalapack.c
        ${CMAKE_SOURCE_DIR}/src/matrixlapack_elpa.c)
    target_include_directories(elpa_statepanel_check PRIVATE
        ${CMAKE_SOURCE_DIR}/src/include ${ELPA_INCLUDE_DIRS})
    target_compile_options(elpa_statepanel_check PRIVATE ${ELPA_COMPILE_OPTIONS})
    target_compile_definitions(elpa_statepanel_check PRIVATE -D_ELPA -D_SCALAPACK -DMPI)
    # _ec_sc_libs is the normalized SCALAPACK_LIBRARIES list computed above
    # for elpa_eigen_check; reuse it here for the same reason.
    target_link_libraries(elpa_statepanel_check ${ELPA_LIBRARIES}
        ${_ec_sc_libs} ${LAPACK_LIBRARIES} ${MPI_C_LIBRARIES}
        ${MPI_Fortran_LIBRARIES} ${HPHI_FORTRAN_RUNTIME_LIBRARIES} m)
    # Same run_with_mpi_precheck.sh + SKIP_RETURN_CODE 77 convention as
    # elpa_redist_check; RedistBlockCyclicToStatePanel's default check
    # (N=97) also builds its BLACS grid from MPI_Dims_create and works
    # with any process count, so "min:2" is enough.
    add_test(
        NAME elpa_statepanel_check
        COMMAND /bin/sh ${CMAKE_SOURCE_DIR}/test/run_with_mpi_precheck.sh
                min:2 ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                ${MPIEXEC_PREFLAGS} $<TARGET_FILE:elpa_statepanel_check>
                ${MPIEXEC_POSTFLAGS}
    )
    set_tests_properties(
        elpa_statepanel_check
        PROPERTIES SKIP_RETURN_CODE 77 LABELS "fulldiag;elpa;unit"
    )

    # Dedicated small-N edge case: reuse the same source (NDIM
    # overridable at compile time) with N=2 and exactly 3 ranks, so the
    # frozen ownership formula (NC=ceil(2/3)=1) guarantees rank 2 owns
    # zero states -- exercising the "zero-owner ranks participate with a
    # 1-element dummy buffer + valid descriptor" path that is easy to
    # miss with only the np-vs-N=97 combinations exercised above.
    add_executable(elpa_statepanel_zero_owner unit/elpa_statepanel_check.c
        ${CMAKE_SOURCE_DIR}/src/matrixscalapack.c
        ${CMAKE_SOURCE_DIR}/src/matrixlapack_elpa.c)
    target_include_directories(elpa_statepanel_zero_owner PRIVATE
        ${CMAKE_SOURCE_DIR}/src/include ${ELPA_INCLUDE_DIRS})
    target_compile_options(elpa_statepanel_zero_owner PRIVATE ${ELPA_COMPILE_OPTIONS})
    target_compile_definitions(elpa_statepanel_zero_owner PRIVATE -D_ELPA -D_SCALAPACK -DMPI -DNDIM=2)
    target_link_libraries(elpa_statepanel_zero_owner ${ELPA_LIBRARIES}
        ${_ec_sc_libs} ${LAPACK_LIBRARIES} ${MPI_C_LIBRARIES}
        ${MPI_Fortran_LIBRARIES} ${HPHI_FORTRAN_RUNTIME_LIBRARIES} m)
    add_test(
        NAME elpa_statepanel_zero_owner
        COMMAND /bin/sh ${CMAKE_SOURCE_DIR}/test/run_with_mpi_precheck.sh
                exact:3 ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3
                ${MPIEXEC_PREFLAGS} $<TARGET_FILE:elpa_statepanel_zero_owner>
                ${MPIEXEC_POSTFLAGS}
    )
    set_tests_properties(
        elpa_statepanel_zero_owner
        PROPERTIES SKIP_RETURN_CODE 77 LABELS "fulldiag;elpa;unit"
    )

    # ExpecMode 0/1/(2) equivalence (phase 3a Task 7): same FullDiag problem,
    # ExpecMode 0 (serial) vs ExpecMode 1 (state-task-parallel) vs a single
    # ExpecMode 2 sanity case (3a: 2 downgrades to 1, so it must reproduce
    # ExpecMode 1's own result and print the downgrade INFO), compared on
    # both zvo_phys_* (S^2/Sz included) and every Green aggregate file
    # (per-state and aggregated kinds). Requires Solver 1 (ScaLAPACK) or 3
    # (ELPA) for ExpecMode to be accepted at all (see readdef.c's
    # cErrExpecMode gate), so -- like fulldiag_elpa_hubbard_chain above --
    # this is genuinely registration-only on a non-ELPA/non-ScaLAPACK build:
    # there is no way to exercise it without one of those solvers compiled
    # in, regardless of MPI availability.
    # Registered as two SEPARATE named ctest cases (not one min:2 case)
    # because np=3 specifically exercises non-divisible state-panel
    # ownership (script takes its process count from ${MPIRUN} via the
    # precheck wrapper; it does not launch a second np internally). Both
    # names point at the SAME script (test/fulldiag_expecmode_equiv.sh),
    # so add_test() is used directly here rather than the
    # add_hphi_mpi_test() macro (which assumes the ctest NAME and the
    # script's base filename match).
    add_test(
        NAME fulldiag_expecmode_equiv_np2
        COMMAND /bin/sh ${CMAKE_SOURCE_DIR}/test/run_with_mpi_precheck.sh
                exact:2 ${CMAKE_SOURCE_DIR}/test/fulldiag_expecmode_equiv.sh
    )
    add_test(
        NAME fulldiag_expecmode_equiv_np3
        COMMAND /bin/sh ${CMAKE_SOURCE_DIR}/test/run_with_mpi_precheck.sh
                exact:3 ${CMAKE_SOURCE_DIR}/test/fulldiag_expecmode_equiv.sh
    )
    set_tests_properties(
        fulldiag_expecmode_equiv_np2 fulldiag_expecmode_equiv_np3
        PROPERTIES SKIP_RETURN_CODE 77 LABELS "fulldiag;elpa;expecmode"
    )
endif(USE_ELPA)

# This case exercises the ScaLAPACK solver and must remain available in a
# USE_SCALAPACK=ON, USE_ELPA=OFF build as well as in an ELPA build.
if(USE_SCALAPACK)
    add_test(
        NAME fulldiag_expecmode_scalapack_np2
        COMMAND /bin/sh ${CMAKE_SOURCE_DIR}/test/run_with_mpi_precheck.sh
                exact:2 ${CMAKE_SOURCE_DIR}/test/fulldiag_expecmode_equiv.sh
    )
    set_tests_properties(
        fulldiag_expecmode_scalapack_np2
        PROPERTIES SKIP_RETURN_CODE 77 LABELS "fulldiag;scalapack;expecmode"
    )
    set_property(
        TEST fulldiag_expecmode_scalapack_np2
        APPEND PROPERTY ENVIRONMENT "HPHI_FULLDIAG_SOLVER=1")
endif(USE_SCALAPACK)
