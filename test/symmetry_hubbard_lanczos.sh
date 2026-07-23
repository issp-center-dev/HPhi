#!/bin/sh -e

mkdir -p symmetry_hubbard_lanczos
cd symmetry_hubbard_lanczos

write_calcmod() {
cat > calcmod.def <<EOF
CalcType 0
CalcModel 0
OutputMode 0
CalcEigenVec 0
InitialVecType 0
OutputEigenVec 0
InputEigenVec 0
OutputHam 0
InputHam 0
ReStart 0
CalcSpec 0
EOF
}

write_modpara() {
cat > modpara.def <<EOF
--------------------
Model_Parameters 0
--------------------
--------------------
--------------------
CDataFileHead zvo
CParaFileHead zqp
--------------------
Nsite 4
Nup 1
Ndown 1
Lanczos_max 20
initial_iv -1
exct 1
LanczosEps 12
LanczosTarget 1
LargeValue 50
EOF
}

write_locspn() {
cat > locspn.def <<EOF
================================
NlocalSpin     0
================================
========i_0LocSpn_Sr=Sr_i=======
================================
EOF
}

write_transfer_ring() {
cat > transfer.def <<EOF
================
NTransfer 16
================
========i s j t t_ij======
================
0 0 1 0 1.0 0.0
1 0 0 0 1.0 0.0
0 1 1 1 1.0 0.0
1 1 0 1 1.0 0.0
1 0 2 0 1.0 0.0
2 0 1 0 1.0 0.0
1 1 2 1 1.0 0.0
2 1 1 1 1.0 0.0
2 0 3 0 1.0 0.0
3 0 2 0 1.0 0.0
2 1 3 1 1.0 0.0
3 1 2 1 1.0 0.0
3 0 0 0 1.0 0.0
0 0 3 0 1.0 0.0
3 1 0 1 1.0 0.0
0 1 3 1 1.0 0.0
EOF
}

write_coulombintra() {
cat > coulombintra.def <<EOF
================
NCoulombIntra 4
================
========i U_i ======
================
0 0.5
1 0.5
2 0.5
3 0.5
EOF
}

write_coulombinter() {
cat > coulombinter.def <<EOF
================
NCoulombInter 1
================
========i_j_V ======
================
0 1 0.25
EOF
}

write_k0_transsym() {
cat > qptransidx.def <<EOF
=============================================
NQPTrans          4
=============================================
======== TrIdx_TrWeight_and_TrIdx_i_xi ======
=============================================
0 1.0
1 1.0
2 1.0
3 1.0
0 0 0 1
0 1 1 1
0 2 2 1
0 3 3 1
1 0 1 1
1 1 2 1
1 2 3 1
1 3 0 1
2 0 2 1
2 1 3 1
2 2 0 1
2 3 1 1
3 0 3 1
3 1 0 1
3 2 1 1
3 3 2 1
EOF
}

write_kpi2_transsym() {
cat > qptransidx.def <<EOF
=============================================
NQPTrans          4
=============================================
======== TrIdx_TrWeight_and_TrIdx_i_xi ======
=============================================
0 1.0 0.0
1 0.0 -1.0
2 -1.0 0.0
3 0.0 1.0
0 0 0 1
0 1 1 1
0 2 2 1
0 3 3 1
1 0 1 1
1 1 2 1
1 2 3 1
1 3 0 1
2 0 2 1
2 1 3 1
2 2 0 1
2 3 1 1
3 0 3 1
3 1 0 1
3 2 1 1
3 3 2 1
EOF
}

write_ref_namelist() {
cat > namelist.def <<EOF
CalcMod calcmod.def
ModPara modpara.def
LocSpin locspn.def
Trans transfer.def
CoulombIntra coulombintra.def
EOF
}

write_sym_namelist() {
    with_coulomb="$1"
cat > namelist.def <<EOF
CalcMod calcmod.def
ModPara modpara.def
LocSpin locspn.def
Trans transfer.def
TransSym qptransidx.def
EOF
    if [ "${with_coulomb}" = "yes" ]; then
cat >> namelist.def <<EOF
CoulombIntra coulombintra.def
EOF
    fi
}

assert_energy() {
    expected="$1"
    log="$2"
    energy=`awk '$1 == "Energy" {print $2; exit}' output/zvo_energy.dat`
    test -n "${energy}"
    diff=`awk -v a="${energy}" -v b="${expected}" 'BEGIN{d=a-b; if(d<0)d=-d; printf "%8.6f", d}'`
    if [ "${diff}" != "0.000000" ]; then
        cat "${log}"
        echo "Energy mismatch: got ${energy}, expected ${expected}"
        exit 1
    fi
}

assert_energy_matches_reference() {
    expected="$1"
    log="$2"
    energy=`awk '$1 == "Energy" {print $2; exit}' output/zvo_energy.dat`
    test -n "${energy}"
    diff=`awk -v a="${energy}" -v b="${expected}" 'BEGIN{d=a-b; if(d<0)d=-d; printf "%8.6f", d}'`
    if [ "${diff}" != "0.000000" ]; then
        cat "${log}"
        echo "Energy mismatch: got ${energy}, reference ${expected}"
        exit 1
    fi
}

assert_doublon_matches_reference() {
    expected="$1"
    log="$2"
    doublon=`awk '$1 == "Doublon" {print $2; exit}' output/zvo_energy.dat`
    test -n "${doublon}"
    if ! awk -v a="${doublon}" -v b="${expected}" \
        'BEGIN{d=a-b; if(d<0)d=-d; exit(d <= 1.0e-6 ? 0 : 1)}'; then
        cat "${log}"
        echo "Doublon mismatch: got ${doublon}, reference ${expected}"
        exit 1
    fi
}

expect_failure() {
    pattern="$1"
    log="$2"
    shift 2
    if "$@" > "${log}" 2>&1; then
        cat "${log}"
        exit 1
    fi
    if ! grep -q "${pattern}" "${log}"; then
        cat "${log}"
        echo "Expected pattern not found: ${pattern}"
        exit 1
    fi
}

assert_symmetry_log() {
    expected_dim="$1"
    log="$2"
    if ! grep -q "Symmetry basis: raw_dim=16 sector_dim=${expected_dim} group_order=4" "${log}"; then
        cat "${log}"
        echo "Expected Hubbard symmetry sector_dim=${expected_dim} was not found"
        exit 1
    fi
    if grep -q "MPI site separation summary" "${log}"; then
        cat "${log}"
        echo "TransSym Hubbard MPI path unexpectedly used site decomposition."
        exit 1
    fi
}

assert_rank_stats() {
    expected_dim="$1"
    expected_ranks="$2"
    log="$3"
    expected_digest="${4:-}"
    expected_reference="${5:-0}"
    stats=output/CalcTimerRankStats.dat
    if [ ! -f output/CalcTimer.dat ]; then
        return
    fi
    if [ ! -f "${stats}" ]; then
        cat "${log}"
        echo "Missing ${stats}"
        exit 1
    fi
    if ! awk -v expected_dim="${expected_dim}" -v expected_ranks="${expected_ranks}" \
        -v expected_digest="${expected_digest}" \
        -v expected_reference="${expected_reference}" '
        function abs(x) { return x < 0 ? -x : x }
        function value(field, parts) {
            split(field, parts, "=")
            return parts[2]
        }
        /^format=/ {
            header_version = value($2)
            header_ranks = value($3)
            header_basis_layout = value($4)
            header_matvec_mode = value($5)
            header_vector_exchange = value($6)
            next
        }
        $1 == "timer" {
            timer_count++
            id = value($2)
            ranks = value($3)
            min = value($4)
            max = value($5)
            mean = value($6)
            if (ranks != expected_ranks || min > mean || mean > max) bad = 1
            if (expected_ranks == 1 &&
                (abs(min - max) > 1.0e-15 || abs(min - mean) > 1.0e-15)) bad = 1
            timer_mean[id] = mean
            next
        }
        $1 == "work" {
            work_count++
            key = value($2)
            ranks = value($3)
            min = value($4)
            max = value($5)
            mean = value($6)
            if (ranks != expected_ranks || min > mean || mean > max) bad = 1
            work_min[key] = min
            work_max[key] = max
            work_mean[key] = mean
            next
        }
        $1 == "metric" {
            metric_count++
            key = value($2)
            ranks = value($3)
            min = value($4)
            max = value($5)
            mean = value($6)
            if (ranks != expected_ranks || min > mean || mean > max) bad = 1
            metric_min[key] = min
            metric_max[key] = max
            metric_mean[key] = mean
            next
        }
        $1 == "basis_digest" {
            digest_count++
            ranks = value($3)
            digest_min = value($4)
            digest_max = value($5)
            if (ranks != expected_ranks || digest_min != digest_max) bad = 1
            if (expected_digest != "" && digest_min != expected_digest) bad = 1
            next
        }
        $1 == "halo_schedule_digest" {
            schedule_digest_count++
            ranks = value($3)
            schedule_digest_xor = value($4)
            schedule_digest_sum = value($5)
            if (ranks != expected_ranks ||
                schedule_digest_xor == "" || schedule_digest_sum == "") bad = 1
            next
        }
        END {
            if (header_version != 3 || header_ranks != expected_ranks ||
                header_basis_layout != "replicated" ||
                header_matvec_mode != "plan" ||
                header_vector_exchange != "allgather" ||
                timer_count != 22 || work_count != 36 ||
                metric_count != 13 || digest_count != 1 ||
                schedule_digest_count != 1) bad = 1
            if (abs(work_mean["basis_raw_states"] * expected_ranks - 16) > 1.0e-12) bad = 1
            if (abs(work_mean["basis_representative_candidates"] * expected_ranks - 4) > 1.0e-12) bad = 1
            if (abs(work_mean["basis_compatible_survivors"] * expected_ranks - expected_dim) > 1.0e-12) bad = 1
            if (work_min["basis_transform_calls"] <= 0) bad = 1
            if (abs(work_mean["basis_orbit_metadata_calls"] * expected_ranks - 4) > 1.0e-12) bad = 1
            if (work_min["basis_thread_count"] < 1) bad = 1
            if (abs(work_mean["plan_local_rows"] * expected_ranks - expected_dim) > 1.0e-12) bad = 1
            if (abs(work_mean["plan_local_column_nnz"] + work_mean["plan_remote_column_nnz"] - work_mean["plan_local_nnz"]) > 1.0e-12) bad = 1
            if (work_max["halo_ghost_count"] > work_max["plan_remote_column_nnz"]) bad = 1
            if (work_max["halo_incoming_peer_count"] >= expected_ranks ||
                work_max["halo_outgoing_peer_count"] >= expected_ranks) bad = 1
            if (work_min["column_slot_width"] != 32 ||
                work_max["column_slot_width"] != 32) bad = 1
            if (work_min["halo_schedule_ready"] != 1 ||
                work_max["halo_schedule_ready"] != 1) bad = 1
            if (work_min["halo_schedule_bytes"] <= 0) bad = 1
            if (abs(work_mean["halo_runtime_buffer_bytes"] - 16 * (work_mean["halo_ghost_count"] + work_mean["halo_send_value_count"])) > 1.0e-12) bad = 1
            if (work_min["halo_reference_enabled"] != expected_reference ||
                work_max["halo_reference_enabled"] != expected_reference) bad = 1
            if (work_min["symmetry_matvec_calls"] <= 0) bad = 1
            if (work_min["prdct_allreduce_calls"] <= 0) bad = 1
            if (abs(work_mean["symmetry_matvec_calls"] - work_mean["prdct_allreduce_calls"]) > 1.0e-12) bad = 1
            if (expected_reference == 1) {
                if (work_min["halo_reference_exchange_calls"] <= 0) bad = 1
                if (abs(work_mean["halo_reference_exchange_calls"] - work_mean["symmetry_matvec_calls"]) > 1.0e-12) bad = 1
            } else if (work_max["halo_reference_exchange_calls"] != 0) {
                bad = 1
            }
            if (abs(work_mean["input_allgather_payload_bytes_per_call"] - 16 * work_mean["input_allgather_nonlocal_values_per_call"]) > 1.0e-12) bad = 1
            if (expected_ranks > 1) {
                if (abs(work_mean["halo_send_value_count"] - work_mean["halo_ghost_count"]) > 1.0e-12) bad = 1
                if (work_min["input_allgather_calls"] <= 0) bad = 1
                if (abs(work_mean["input_allgather_calls"] - work_mean["symmetry_matvec_calls"]) > 1.0e-12) bad = 1
            } else {
                if (work_max["plan_remote_column_nnz"] != 0 ||
                    work_max["halo_ghost_count"] != 0 ||
                    work_max["input_allgather_calls"] != 0) bad = 1
            }
            if (metric_min["plan_remote_column_nnz_ratio"] < 0 ||
                metric_max["plan_remote_column_nnz_ratio"] > 1 ||
                metric_min["halo_ghost_global_ratio"] < 0 ||
                metric_max["halo_ghost_global_ratio"] > 1 ||
                metric_min["halo_ghost_nonlocal_ratio"] < 0 ||
                metric_max["halo_ghost_nonlocal_ratio"] > 1 ||
                metric_min["symmetry_matvec_seconds_per_call"] < 0 ||
                metric_min["input_allgather_seconds_per_call"] < 0 ||
                metric_min["input_allgather_effective_bandwidth_Bps"] < 0 ||
                metric_min["plan_apply_seconds_per_call"] < 0 ||
                metric_min["prdct_allreduce_seconds_per_call"] < 0 ||
                metric_min["halo_reference_pack_seconds_per_call"] < 0 ||
                metric_min["halo_reference_exchange_seconds_per_call"] < 0 ||
                metric_min["halo_reference_validation_seconds_per_call"] < 0) bad = 1
            exit bad
        }
    ' "${stats}"; then
        cat "${log}"
        cat "${stats}"
        echo "Invalid rank-aware symmetry setup statistics"
        exit 1
    fi
}

run_mpi_symmetry_case() {
    label="$1"
    expected_energy="$2"
    expected_dim="$3"
    expected_doublon="${4:-}"
    expected_ranks="$5"
    expected_digest="$6"
    log_file="hubbard_${label}_mpi.log"
    rm -rf output
    if ! env HPHI_SYMMETRY_HALO_REFERENCE=1 \
        ${MPIRUN} ../../src/HPhi -e namelist.def > "${log_file}" 2>&1; then
        cat "${log_file}"
        exit 1
    fi
    assert_energy_matches_reference "${expected_energy}" "${log_file}"
    if [ -n "${expected_doublon}" ]; then
        assert_doublon_matches_reference "${expected_doublon}" "${log_file}"
    fi
    assert_symmetry_log "${expected_dim}" "${log_file}"
    assert_rank_stats "${expected_dim}" "${expected_ranks}" "${log_file}" \
        "${expected_digest}" 1
}

run_mpi_if_available() {
    label="$1"
    expected_energy="$2"
    expected_dim="$3"
    expected_doublon="${4:-}"
    if [ -n "${MPIRUN}" ]; then
        expected_digest=`awk '$1 == "basis_digest" {
            split($4, parts, "="); print parts[2]; exit
        }' output/CalcTimerRankStats.dat`
        if [ -z "${expected_digest}" ]; then
            echo "Missing serial symmetry basis digest"
            exit 1
        fi
        MPI_NP=`printf "%s\n" "${MPIRUN}" | awk '{for(i=1;i<=NF;i++){if($i=="-np"||$i=="-n"){print $(i+1); exit}}}'`
        if printf "%s\n" "${MPI_NP}" | grep -Eq "^[0-9]+$" && [ "${MPI_NP}" -gt 1 ]; then
            run_mpi_symmetry_case "$label" "$expected_energy" "$expected_dim" \
                "$expected_doublon" "$MPI_NP" "$expected_digest"
        fi
    fi
}

write_calcmod
write_modpara
write_locspn
write_transfer_ring
write_coulombintra
write_ref_namelist

../../src/HPhi -e namelist.def > hubbard_ref.log 2>&1
ref_energy=`awk '$1 == "Energy" {print $2; exit}' output/zvo_energy.dat`
ref_doublon=`awk '$1 == "Doublon" {print $2; exit}' output/zvo_energy.dat`
test -n "${ref_energy}"
test -n "${ref_doublon}"
rm -rf output

write_k0_transsym
write_sym_namelist yes
../../src/HPhi -e namelist.def > hubbard_k0.log 2>&1
assert_energy_matches_reference "${ref_energy}" hubbard_k0.log
assert_doublon_matches_reference "${ref_doublon}" hubbard_k0.log
grep -q "Symmetry basis: raw_dim=16 sector_dim=4 group_order=4" hubbard_k0.log
if grep -q "MPI site separation summary" hubbard_k0.log; then
    cat hubbard_k0.log
    echo "TransSym Hubbard serial path unexpectedly used site decomposition."
    exit 1
fi
assert_rank_stats 4 1 hubbard_k0.log
run_mpi_if_available k0 "${ref_energy}" 4 "${ref_doublon}"
expect_failure "HPHI_SYMMETRY_HALO_REFERENCE must be" \
    invalid_halo_reference.log env HPHI_SYMMETRY_HALO_REFERENCE=invalid \
    ../../src/HPhi -e namelist.def

rm -rf output
write_kpi2_transsym
write_sym_namelist no
../../src/HPhi -e namelist.def > hubbard_kpi2.log 2>&1
assert_energy "-2.0" hubbard_kpi2.log
grep -q "Symmetry basis: raw_dim=16 sector_dim=4 group_order=4" hubbard_kpi2.log
if grep -q "MPI site separation summary" hubbard_kpi2.log; then
    cat hubbard_kpi2.log
    echo "TransSym Hubbard serial path unexpectedly used site decomposition."
    exit 1
fi
assert_rank_stats 4 1 hubbard_kpi2.log
run_mpi_if_available kpi2 "-2.0" 4

rm -rf output
write_calcmod
write_k0_transsym
write_sym_namelist yes
perl -0pi -e 's/CalcType 0/CalcType 3/' calcmod.def
env HPHI_SYMMETRY_HALO_REFERENCE=1 \
    ../../src/HPhi -e namelist.def > hubbard_k0_cg.log 2>&1
assert_energy_matches_reference "${ref_energy}" hubbard_k0_cg.log
assert_doublon_matches_reference "${ref_doublon}" hubbard_k0_cg.log
assert_symmetry_log 4 hubbard_k0_cg.log
assert_rank_stats 4 1 hubbard_k0_cg.log "" 1
run_mpi_if_available k0_cg "${ref_energy}" 4 "${ref_doublon}"
write_calcmod

rm -rf output
write_k0_transsym
write_sym_namelist yes
perl -0pi -e 's/2 0 3 0 1.0 0.0/2 0 3 0 0.5 0.0/' transfer.def
perl -0pi -e 's/3 0 2 0 1.0 0.0/3 0 2 0 0.5 0.0/' transfer.def
expect_failure "Hubbard Transfer invariance failed" \
    noninvariant_transfer.log ../../src/HPhi -e namelist.def

write_transfer_ring
perl -0pi -e 's/2 0.5/2 0.25/' coulombintra.def
expect_failure "Hubbard CoulombIntra invariance failed" \
    noninvariant_coulombintra.log ../../src/HPhi -e namelist.def

write_coulombintra
write_coulombinter
cat >> namelist.def <<EOF
CoulombInter coulombinter.def
EOF
expect_failure "Hubbard symmetry basis supports Transfer and CoulombIntra terms only" \
    unsupported_term.log ../../src/HPhi -e namelist.def

if [ -n "${MPIRUN}" ]; then
    MPI_NP=`printf "%s\n" "${MPIRUN}" | awk '{for(i=1;i<=NF;i++){if($i=="-np"||$i=="-n"){print $(i+1); exit}}}'`
    if printf "%s\n" "${MPI_NP}" | grep -Eq "^[0-9]+$" && [ "${MPI_NP}" -gt 1 ]; then
        expect_failure "Hubbard symmetry basis supports Transfer and CoulombIntra terms only" \
            unsupported_term_mpi.log ${MPIRUN} ../../src/HPhi -e namelist.def
    fi
fi
