#!/bin/sh -e
set -e

# Reflection symmetry of a SpinlessFermion ring through the expert-mode
# TransSym file.  A reflection permutes the occupied orbitals, so the fermion
# sign of the permutation matters: with 3 fermions on 6 sites, the 4 states
# that the site-centred mirror R(i) = -i maps onto themselves all contain one
# exchanged pair and therefore have R = -1.  The two parity sectors then have
# dimensions 8 (R = +1) and 12 (R = -1); a bosonic counting would give 12 and 8.
# The bond-centred mirror R(i) = 5 - i fixes no state (10 and 10).
#
# The parity sectors together must reproduce the complete raw spectrum
# (20 levels) and their dimensions must add up to 20.  The raw reference is a
# CG run with exct = 20, because FullDiag does not support SpinlessFermion.
# The sector dimensions were also confirmed independently of HPhi from
# signed permutation traces.

mkdir -p symmetry_spinless_fermion_point_group
cd symmetry_spinless_fermion_point_group

L=6
NE=3
RAW_DIM=20

run_hphi() {
    log="$1"
    shift
    "$@" > "${log}" 2>&1 || { cat "${log}"; exit 1; }
}

cat > calcmod.def <<EOF
CalcType 3
CalcModel 7
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
Nsite ${L}
Ncond ${NE}
Lanczos_max 400
initial_iv 1
exct $1
LanczosEps 12
LanczosTarget 2
LargeValue 50
PreCG 0
EOF
}

cat > locspn.def <<EOF
================================
NlocalSpin     0
================================
========i_0LocSpn_Sr=Sr_i=======
================================
EOF

# Nearest-neighbour hopping t = 1 and density-density V = 0.5 on the ring.
awk -v L="${L}" 'BEGIN {
    print "================"; printf "NTransfer %d\n", 2 * L
    print "================"; print "========i s j t t_ij======"; print "================"
    for (i = 0; i < L; i++) {
        j = (i + 1) % L
        printf "%d 0 %d 0 1.0 0.0\n", i, j
        printf "%d 0 %d 0 1.0 0.0\n", j, i
    }
}' > transfer.def
awk -v L="${L}" 'BEGIN {
    print "================"; printf "NCoulombInter %d\n", L
    print "================"; print "========i_j_V ======"; print "================"
    for (i = 0; i < L; i++) printf "%d %d 0.5\n", i, (i + 1) % L
}' > coulombinter.def

cat > namelist_ref.def <<EOF
CalcMod calcmod.def
ModPara modpara.def
LocSpin locspn.def
Trans transfer.def
CoulombInter coulombinter.def
EOF

cat > namelist.def <<EOF
CalcMod calcmod.def
ModPara modpara.def
LocSpin locspn.def
Trans transfer.def
CoulombInter coulombinter.def
TransSym qptransidx.def
EOF

# write_mirror <chi_R> <c>: group {e, R}, R(i) = (c - i) mod L, chi(R) = chi_R.
write_mirror() {
    awk -v L="${L}" -v sr="$1" -v c="$2" 'BEGIN {
        print "============================================="
        print "NQPTrans          2"
        print "============================================="
        print "======== TrIdx_TrWeight_and_TrIdx_i_xi ======"
        print "============================================="
        print "0 1"
        printf "1 %d\n", sr
        for (i = 0; i < L; i++) printf "0 %d %d 1\n", i, i
        for (i = 0; i < L; i++) printf "1 %d %d 1\n", i, (c - i + L) % L
    }' > qptransidx.def
}

sorted_energies() {
    awk '$1 == "Energy" {print $2}' output/zvo_energy.dat | sort -g
}

run_sector() {
    label="$1"
    dim="$2"
    write_modpara "${dim}"
    rm -rf output
    run_hphi "${label}.log" ${MPIRUN} ../../src/HPhi -e namelist.def
    grep -q "Symmetry basis: raw_dim=${RAW_DIM} sector_dim=${dim} group_order=2" "${label}.log"
    sorted_energies > "${label}.energies"
    test "`wc -l < "${label}.energies" | tr -d ' '`" = "${dim}"
}

assert_spectrum() {
    count_a=`wc -l < "$2" | tr -d ' '`
    count_b=`wc -l < "$3" | tr -d ' '`
    if [ "${count_a}" != "${count_b}" ]; then
        echo "[$1] level count mismatch: ${count_a} vs ${count_b}" >&2
        exit 1
    fi
    max_diff=`paste "$2" "$3" | awk 'BEGIN{m=0} {d=$1-$2; if(d<0)d=-d; if(d>m)m=d} END{printf "%.3e", m}'`
    if ! awk -v d="${max_diff}" -v tol="$4" 'BEGIN{exit !(d <= tol)}'; then
        echo "[$1] spectra differ: max |diff| = ${max_diff} > $4" >&2
        exit 1
    fi
    echo "[$1] ${count_a} levels agree (max |diff| = ${max_diff})"
}

# Raw reference: all 20 levels (serial, no TransSym).
write_modpara ${RAW_DIM}
rm -rf output
run_hphi raw.log ../../src/HPhi -e namelist_ref.def
sorted_energies > raw.energies
test "`wc -l < raw.energies | tr -d ' '`" = "${RAW_DIM}"

# check_mirror <name> <offset c> <dim of chi_R=+1> <dim of chi_R=-1>
check_mirror() {
    name="$1"
    write_mirror 1 "$2"
    run_sector "mirror_${name}_even" "$3"
    write_mirror -1 "$2"
    run_sector "mirror_${name}_odd" "$4"
    test $(($3 + $4)) -eq ${RAW_DIM}
    sort -g "mirror_${name}_even.energies" "mirror_${name}_odd.energies" > "mirror_${name}.union"
    assert_spectrum "spinless mirror ${name}" "mirror_${name}.union" raw.energies 1.0e-8
}

check_mirror site_centred 0 8 12
check_mirror bond_centred 5 10 10

exit 0
