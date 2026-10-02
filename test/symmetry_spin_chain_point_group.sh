#!/bin/sh -e
set -e

# Non-translation (point-group) symmetries through the expert-mode TransSym
# file.  TransSym accepts any permutation group of the sites that has a
# one-dimensional character, not only translations.  This test checks
# reflections and the dihedral group (translation x reflection) of periodic
# Heisenberg rings (Spin-1/2, fixed 2Sz = 0).
#
#   (A) Mirror group {e, R} of a 6-site ring, for a site-centred mirror
#       R(i) = -i and a bond-centred mirror R(i) = 5 - i.  The two parity
#       sectors together must reproduce the complete raw spectrum (20 levels,
#       FullDiag) and their dimensions must add up to 20.
#   (B) Dihedral group D8 (order 16) of an 8-site ring.  Its one-dimensional
#       irreps (k, parity) with k = 0, pi split each translation-only sector
#       in two: the union of the two parity sectors must reproduce the
#       complete translation-only spectrum of the same k.
#   (C) A Hamiltonian that is not invariant under the mirror (one bond with a
#       different coupling) must be rejected at input.
#
# The expected sector dimensions follow from character counting, e.g.
# (20 +- 4)/2 for the site-centred mirror, which fixes 4 of the 20 states.
# They were confirmed independently of HPhi from permutation traces, and
# every sector is solved completely (exct = sector dimension).

mkdir -p symmetry_spin_chain_point_group
cd symmetry_spin_chain_point_group

run_hphi() {
    log="$1"
    shift
    "$@" > "${log}" 2>&1 || { cat "${log}"; exit 1; }
}

write_calcmod() {
cat > calcmod.def <<EOF
CalcType $1
CalcModel 1
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
Nsite $1
2Sz 0
Lanczos_max 400
initial_iv 1
exct $2
LanczosEps 12
LanczosTarget 2
LargeValue 50
PreCG 0
EOF
}

# write_ring <L> <bond0_coupling>: Heisenberg ring, all bonds J = 1 except
# the bond (0,1), whose Exchange coupling is <bond0_coupling>.
write_ring() {
    awk -v L="$1" -v j0="$2" 'BEGIN {
        print "================"; printf "NlocalSpin %d\n", L
        print "================"; print "========i_1LocSpn ======"; print "================"
        for (i = 0; i < L; i++) printf "%d 1\n", i
    }' > locspn.def
    awk -v L="$1" -v j0="$2" 'BEGIN {
        print "================"; printf "NExchange %d\n", L
        print "================"; print "========i_j_J ======"; print "================"
        for (i = 0; i < L; i++) printf "%d %d %s\n", i, (i + 1) % L, (i == 0 ? j0 : "1.0")
    }' > exchange.def
    awk -v L="$1" 'BEGIN {
        print "================"; printf "NIsing %d\n", L
        print "================"; print "========i_j_J ======"; print "================"
        for (i = 0; i < L; i++) printf "%d %d 1.0\n", i, (i + 1) % L
    }' > ising.def
}

cat > namelist_ref.def <<EOF
CalcMod calcmod.def
ModPara modpara.def
LocSpin locspn.def
Exchange exchange.def
Ising ising.def
EOF

cat > namelist.def <<EOF
CalcMod calcmod.def
ModPara modpara.def
LocSpin locspn.def
Exchange exchange.def
Ising ising.def
TransSym qptransidx.def
EOF

# write_group <kind> <L> <chi_T> <chi_R> <c>: qptransidx.def of
#   mirror      {e, R},       R(i) = (c - i) mod L, character chi_R^b
#   translation {T^a},        T(i) = i + 1,         character chi_T^a
#   dihedral    {T^a R^b},    g(i) = a + (-1)^b i,  character chi_T^a chi_R^b
# Operation n = b * (number of translations) + a.  chi_T = +-1 (L even).
write_group() {
    awk -v kind="$1" -v L="$2" -v st="$3" -v sr="$4" -v c="$5" 'BEGIN {
        if (kind == "mirror") { na = 1; nb = 2 }
        else if (kind == "translation") { na = L; nb = 1 }
        else { na = L; nb = 2 }
        print "============================================="
        printf "NQPTrans          %d\n", na * nb
        print "============================================="
        print "======== TrIdx_TrWeight_and_TrIdx_i_xi ======"
        print "============================================="
        for (b = 0; b < nb; b++) for (a = 0; a < na; a++) {
            ch = 1
            for (k = 0; k < a; k++) ch *= st
            for (k = 0; k < b; k++) ch *= sr
            printf "%d %d\n", b * na + a, ch
        }
        for (b = 0; b < nb; b++) for (a = 0; a < na; a++)
            for (i = 0; i < L; i++) {
                if (b == 0) t = (i + a) % L
                else if (kind == "mirror") t = (c - i + L) % L
                else t = (a - i + L) % L
                printf "%d %d %d 1\n", b * na + a, i, t
            }
    }' > qptransidx.def
}

sorted_energies() {
    awk '$1 == "Energy" {print $2}' output/zvo_energy.dat | sort -g
}

# run_sector <label> <sector dimension> <group order>: solve the complete
# spectrum (exct = dimension) of the sector defined by qptransidx.def.
run_sector() {
    label="$1"
    dim="$2"
    order="$3"
    write_modpara "${L}" "${dim}"
    rm -rf output
    run_hphi "${label}.log" ${MPIRUN} ../../src/HPhi -e namelist.def
    grep -q "Symmetry basis: raw_dim=${RAW_DIM} sector_dim=${dim} group_order=${order}" "${label}.log"
    sorted_energies > "${label}.energies"
    test "`wc -l < "${label}.energies" | tr -d ' '`" = "${dim}"
}

# assert_spectrum <label> <file a> <file b> <tolerance>
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

# ---------------------------------------------------------------------------
# (A) Mirror groups of the 6-site ring against the raw FullDiag spectrum.
L=6
RAW_DIM=20
write_ring ${L} 1.0
write_calcmod 2
write_modpara ${L} 1
rm -rf output
run_hphi raw_l6.log ../../src/HPhi -e namelist_ref.def
awk '{print $2}' output/Eigenvalue.dat | sort -g > raw_l6.energies
test "`wc -l < raw_l6.energies | tr -d ' '`" = "${RAW_DIM}"

write_calcmod 3

# check_mirror <name> <offset c> <dim of chi_R=+1> <dim of chi_R=-1>
check_mirror() {
    name="$1"
    write_group mirror ${L} 1 1 "$2"
    run_sector "mirror_${name}_even" "$3" 2
    write_group mirror ${L} 1 -1 "$2"
    run_sector "mirror_${name}_odd" "$4" 2
    test $(($3 + $4)) -eq ${RAW_DIM}
    sort -g "mirror_${name}_even.energies" "mirror_${name}_odd.energies" > "mirror_${name}.union"
    assert_spectrum "mirror ${name}" "mirror_${name}.union" raw_l6.energies 1.0e-8
}

check_mirror site_centred 0 12 8
check_mirror bond_centred 5 10 10

# ---------------------------------------------------------------------------
# (B) Dihedral group D8 of the 8-site ring against translation-only sectors.
L=8
RAW_DIM=70
write_ring ${L} 1.0

# check_dihedral <name> <chi_T> <dim of translation sector> <dim of chi_R=+1> <dim of chi_R=-1>
check_dihedral() {
    name="$1"
    write_group translation ${L} "$2" 1 0
    run_sector "translation_${name}" "$3" 8
    write_group dihedral ${L} "$2" 1 0
    run_sector "dihedral_${name}_even" "$4" 16
    write_group dihedral ${L} "$2" -1 0
    run_sector "dihedral_${name}_odd" "$5" 16
    test $(($4 + $5)) -eq "$3"
    sort -g "dihedral_${name}_even.energies" "dihedral_${name}_odd.energies" > "dihedral_${name}.union"
    assert_spectrum "D8 ${name}" "dihedral_${name}.union" "translation_${name}.energies" 1.0e-8
}

check_dihedral k0 1 10 8 2
check_dihedral kpi -1 10 5 5

# ---------------------------------------------------------------------------
# (C) A Hamiltonian that the mirror does not leave invariant must be rejected.
L=6
RAW_DIM=20
write_ring ${L} 2.0
write_group mirror ${L} 1 1 0
write_modpara ${L} 1
rm -rf output
if ${MPIRUN} ../../src/HPhi -e namelist.def > not_invariant.log 2>&1; then
    cat not_invariant.log
    echo "A mirror-asymmetric Hamiltonian was accepted." >&2
    exit 1
fi
grep -q "TransSym Hamiltonian invariance failed for Exchange term" not_invariant.log

exit 0
