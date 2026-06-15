#!/bin/sh -e
# SpinlessFermion diagonal InterAll (density-density n_i n_j) must enter the
# Hamiltonian and give the same energy as the equivalent CoulombInter.
#
# Prerelease finding H-3: SetDiagonalInterAll had no SpinlessFermion /
# SpinlessFermionGC case in any of its three site-placement branches, so a
# diagonal InterAll term fell through to `default: return -1`. diagonalcalc()
# ignored that return value, so the term was silently dropped from H (only a
# non-fatal "CalcModel 7 is incorrect" line was printed) and the energy was
# that of the non-interacting system. CoulombInter handled the same density-
# density term correctly, so the two diverged.
#
# $1 = CMAKE_SOURCE_DIR (add_hphi_test_with_srcdir)
SRCDIR="$1"

mkdir -p spinless_interall_diagonal
cd spinless_interall_diagonal

# Generate the SpinlessFermion expert-mode definition files (L=8 chain, 3
# fermions; non-degenerate ground state) with no interaction.
python3 "${SRCDIR}/test/testSpinlessCalc.py" -p "../../src/HPhi" -m "SpinlessFermion" \
  -s 8 -n 3 -V 0.0 > log_generate.txt 2>&1
cp namelist.def namelist.base

run_energy() {  # echoes the ground-state energy
  rm -rf output
  ../../src/HPhi -e namelist.def > "$1" 2>&1
  awk '/Energy/{print $2}' output/zvo_energy.dat
}

# (a) diagonal InterAll: n_0 n_1, n_1 n_2, n_3 n_5 with distinct couplings.
cp namelist.base namelist.def
cat > interall_diag.def <<EOF
=====
NInterAll 3
=====
=====
=====
0 0 0 0 1 0 1 0 2.0 0.0
1 0 1 0 2 0 2 0 1.5 0.0
3 0 3 0 5 0 5 0 1.0 0.0
EOF
printf '  InterAll interall_diag.def\n' >> namelist.def
e_interall=$(run_energy log_interall.txt)

# (b) the equivalent CoulombInter terms.
cp namelist.base namelist.def
cat > coulombinter_eq.def <<EOF
=====
NCoulombInter 3
=====
=====
=====
0 1 2.0
1 2 1.5
3 5 1.0
EOF
printf '  CoulombInter coulombinter_eq.def\n' >> namelist.def
e_coulomb=$(run_energy log_coulomb.txt)

echo "[spinless-interall-diag] InterAll energy = ${e_interall} ; CoulombInter energy = ${e_coulomb}"
d=$(awk -v a="${e_interall}" -v b="${e_coulomb}" 'BEGIN{x=a-b; if(x<0)x=-x; printf "%.3e", x}')
echo "[spinless-interall-diag] |InterAll - CoulombInter| = ${d}"
awk -v d="${d}" 'BEGIN{ exit (d < 1e-9) ? 0 : 1 }' || {
  echo "[spinless-interall-diag] MISMATCH (diagonal InterAll dropped from H for SpinlessFermion)"; exit 1; }
echo "SpinlessFermion diagonal InterAll == CoulombInter within tolerance."
