#!/bin/sh
# FullDiag with the ScaLAPACK solver on a Hilbert space of dimension 1.
# Some ScaLAPACK implementations (Intel MKL) return the 1x1 eigenvector of
# pzheev as 1+i, whose norm is sqrt(2); every expectation value then came out
# twice as large. Each case has a single state with E = -4 * 0.3 = -1.2, and
# must give <H> = -1.2, <N> = Ncond and n_i = 1 on every orbital.
#  - spinless: SpinlessFermion, 4 sites, 4 particles;
#  - hubbard:  Hubbard, 2 sites, 4 electrons (2Sz = 0).

set -e

HPHI="$(cd "$(pwd)/.." && pwd)/src/HPhi"
mkdir -p fulldiag_scalapack_dim1
cd fulldiag_scalapack_dim1

write_case() {  # $1 = dir, $2 = CalcModel, $3 = Nsite, $4 = modpara lines, $5 = spins
  mkdir -p "$1"
  n=$3
  {
    printf "====\nNlocalSpin 0\n====\n====\n====\n"
    i=0; while [ $i -lt $n ]; do echo "$i 0"; i=$((i + 1)); done
  } > "$1/locspn.def"
  rows=""; green=""; cnt=0
  for s in $5; do
    i=0
    while [ $i -lt $n ]; do
      rows="$rows$i $s $i $s 0.3 0.0
"
      green="$green$i $s $i $s
"
      cnt=$((cnt + 1)); i=$((i + 1))
    done
  done
  printf "====\nNTransfer %d\n====\n====\n====\n%s" $cnt "$rows" > "$1/trans.def"
  printf "====\nNCisAjs %d\n====\n====\n====\n%s" $cnt "$green" > "$1/greenone.def"
  cat > "$1/modpara.def" <<EOF
--------------------
Model_Parameters   0
--------------------
HPhi_Cal_Parameters
--------------------
CDataFileHead  zvo
CParaFileHead  zqp
--------------------
Nsite          $n
$4
Lanczos_max    2000
initial_iv     -1
exct           1
LanczosEps     14
LanczosTarget  2
LargeValue     12.0
NumAve         5
ExpecInterval  20
EOF
  cat > "$1/calcmod.def" <<EOF
CalcType        2
CalcModel       $2
ReStart         0
CalcSpec        0
CalcEigenVec    0
InitialVecType  0
InputEigenVec   0
OutputEigenVec  0
Solver          1
EOF
  cat > "$1/namelist.def" <<EOF
ModPara   modpara.def
LocSpin   locspn.def
Trans     trans.def
OneBodyG  greenone.def
CalcMod   calcmod.def
EOF
}

write_case spinless 7 4 "Ncond          4" "0"
write_case hubbard 0 2 "Ncond          4
2Sz            0" "0 1"

for d in spinless hubbard; do
  ( cd "$d" && ${MPIRUN} "${HPHI}" -e namelist.def > log.txt 2>&1 ) || {
    echo "[fulldiag_scalapack_dim1] HPhi failed in $d"; tail -5 "$d/log.txt"; exit 1; }
done

python3 - <<'PYEOF'
import glob, sys
fail = False
for d in ("spinless", "hubbard"):
    phys = glob.glob(d + "/output/zvo_phys_*.dat")
    if len(phys) != 1:
        sys.exit("[fulldiag_scalapack_dim1] %s: no phys file" % d)
    h, n = (float(x) for x in open(phys[0]).read().split("\n")[1].split()[:2])
    dens = [float(l.split()[4]) for l in open(d + "/output/zvo_cisajs_eigen0.dat") if l.strip()]
    print("[fulldiag_scalapack_dim1] %s: <H> = %.12f, <N> = %.12f, n_i = %s" % (d, h, n, dens))
    if abs(h + 1.2) > 1e-10 or abs(n - 4.0) > 1e-10 or any(abs(x - 1.0) > 1e-10 for x in dens):
        print("[fulldiag_scalapack_dim1] %s: expected <H> = -1.2, <N> = 4 and n_i = 1" % d)
        fail = True
sys.exit(1 if fail else 0)
PYEOF
echo "fulldiag_scalapack_dim1: OK"
