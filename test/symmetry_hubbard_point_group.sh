#!/bin/sh -e
set -e

# Reflection symmetry of a Hubbard ring through the expert-mode TransSym file.
# A reflection permutes the occupied orbitals of each spin, so the fermion
# sign of the permutation matters.  With 2 up and 1 down electrons on 4 sites
# the site-centred mirror R(i) = -i fixes 4 of the 24 states: two with the up
# pair {0, 2} (sign +1) and two with the exchanged up pair {1, 3} (sign -1).
# The two parity sectors therefore have dimensions 12 and 12, whereas a
# bosonic counting (all signs +1) would give 14 and 10.  The bond-centred
# mirror R(i) = 3 - i fixes no state and also gives 12 and 12.
#
# The parity sectors together must reproduce the complete raw spectrum
# (24 levels, FullDiag) and their dimensions must add up to 24.  The sector
# dimensions were also confirmed independently of HPhi from signed
# permutation traces.

mkdir -p symmetry_hubbard_point_group
cd symmetry_hubbard_point_group

L=4
RAW_DIM=24

run_hphi() {
    log="$1"
    shift
    "$@" > "${log}" 2>&1 || { cat "${log}"; exit 1; }
}

write_calcmod() {
cat > calcmod.def <<EOF
CalcType $1
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
Nsite ${L}
Nup 2
Ndown 1
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

# Nearest-neighbour hopping t = 1 for both spins and on-site U = 2 on the ring.
awk -v L="${L}" 'BEGIN {
    print "================"; printf "NTransfer %d\n", 4 * L
    print "================"; print "========i s j t t_ij======"; print "================"
    for (i = 0; i < L; i++) {
        j = (i + 1) % L
        for (s = 0; s < 2; s++) {
            printf "%d %d %d %d 1.0 0.0\n", i, s, j, s
            printf "%d %d %d %d 1.0 0.0\n", j, s, i, s
        }
    }
}' > transfer.def
awk -v L="${L}" 'BEGIN {
    print "================"; printf "NCoulombIntra %d\n", L
    print "================"; print "========i U_i ======"; print "================"
    for (i = 0; i < L; i++) printf "%d 2.0\n", i
}' > coulombintra.def

cat > namelist_ref.def <<EOF
CalcMod calcmod.def
ModPara modpara.def
LocSpin locspn.def
Trans transfer.def
CoulombIntra coulombintra.def
EOF

cat > namelist.def <<EOF
CalcMod calcmod.def
ModPara modpara.def
LocSpin locspn.def
Trans transfer.def
CoulombIntra coulombintra.def
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

# Raw reference: all 24 levels by FullDiag (serial, no TransSym).
write_calcmod 2
write_modpara 1
rm -rf output
run_hphi raw.log ../../src/HPhi -e namelist_ref.def
awk '{print $2}' output/Eigenvalue.dat | sort -g > raw.energies
test "`wc -l < raw.energies | tr -d ' '`" = "${RAW_DIM}"

write_calcmod 3

# check_mirror <name> <offset c> <dim of chi_R=+1> <dim of chi_R=-1>
check_mirror() {
    name="$1"
    write_mirror 1 "$2"
    run_sector "mirror_${name}_even" "$3"
    write_mirror -1 "$2"
    run_sector "mirror_${name}_odd" "$4"
    test $(($3 + $4)) -eq ${RAW_DIM}
    sort -g "mirror_${name}_even.energies" "mirror_${name}_odd.energies" > "mirror_${name}.union"
    assert_spectrum "hubbard mirror ${name}" "mirror_${name}.union" raw.energies 1.0e-8
}

check_mirror site_centred 0 12 12
check_mirror bond_centred 3 12 12

exit 0
