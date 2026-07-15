#!/bin/sh
set -e

testname="lanczos_tj_nconserved_nbodyg"
hphi="../../src/HPhi"
tol="0.00000001"

mkdir -p "${testname}"
cd "${testname}"

run_hphi() {
  log="$1"
  shift
  "$@" > "${log}" 2>&1 || { cat "${log}"; exit 1; }
}

cat > namelist.def <<EOF
       ModPara  modpara.def
       LocSpin  locspn.def
         Trans  trans.def
      InterAll  interall.def
      OneBodyG  greenone.def
      TwoBodyG  greentwo.def
 NBodyInterAll  nbodyinterall.def
        NBodyG  nbodyg.def
       CalcMod  calcmod.def
EOF

cat > calcmod.def <<EOF
CalcType        0
CalcModel       9
ReStart         0
CalcSpec        0
CalcEigenVec    0
InitialVecType  0
InputEigenVec   0
OutputEigenVec  0
InputHam        0
OutputHam       0
OutputExVec     0
EOF

cat > modpara.def <<EOF
--------------------
Model_Parameters   0
--------------------
HPhi_Cal_Parameters
--------------------
CDataFileHead  zvo
CParaFileHead  zqp
--------------------
Nsite             4
Ncond             2
Lanczos_max       2000
initial_iv        1
exct              1
LanczosEps        14
LanczosTarget     2
LargeValue        20.0
NumAve            1
ExpecInterval     20
EOF

cat > locspn.def <<EOF
================================
NlocalSpin     0
================================
========i_1LocSpn_0IteElc ======
================================
    0      0
    1      0
    2      0
    3      0
EOF

cat > trans.def <<EOF
========================
NTransfer      12
========================
========i_j_s_tijs======
========================
0 0 1 0 1.0000000000000000 0.0000000000000000
1 0 0 0 1.0000000000000000 0.0000000000000000
0 1 1 1 1.0000000000000000 0.0000000000000000
1 1 0 1 1.0000000000000000 0.0000000000000000
1 0 2 0 1.0000000000000000 0.0000000000000000
2 0 1 0 1.0000000000000000 0.0000000000000000
1 1 2 1 1.0000000000000000 0.0000000000000000
2 1 1 1 1.0000000000000000 0.0000000000000000
2 0 3 0 1.0000000000000000 0.0000000000000000
3 0 2 0 1.0000000000000000 0.0000000000000000
2 1 3 1 1.0000000000000000 0.0000000000000000
3 1 2 1 1.0000000000000000 0.0000000000000000
EOF

cat > interall.def <<EOF
======================
NInterAll      0
======================
========zInterAll=====
======================
EOF

cat > greenone.def <<EOF
========================
NCisAjs 2
========================
========GreenOne========
========================
2 0 2 1
2 1 2 0
EOF

cat > greentwo.def <<EOF
========================
NCisAjsCktAlt 0
========================
========GreenTwo========
========================
EOF

cat > nbodyinterall.def <<EOF
========================
NNBodyInterAll 2
========================
========NBodyInterAll===
========================
1 2 0 2 1 0.1300000000000000 0.0200000000000000
1 2 1 2 0 0.1300000000000000 -0.0200000000000000
EOF

cat > nbodyg.def <<EOF
========================
NNBodyG 2
========================
========NBodyG==========
========================
1 2 0 2 1
1 2 1 2 0
EOF

run_hphi log_lanczos.txt ${MPIRUN} "${hphi}" -e namelist.def
test -f output/zvo_NBodyG.dat || { echo "zvo_NBodyG.dat was not generated"; exit 1; }

awk -v t="${tol}" '
  NR == FNR {
    if ($1 == 2 && $2 == 0 && $3 == 2 && $4 == 1) {
      ref_re = $5;
      ref_im = $6;
      found_ref = 1;
    }
    next;
  }
  $1 == 1 && $2 == 2 && $3 == 0 && $4 == 2 && $5 == 1 {
    dr = $6 - ref_re;
    di = $7 - ref_im;
    if (dr < 0) dr = -dr;
    if (di < 0) di = -di;
    found_nbody = 1;
    ok = (found_ref && dr < t && di < t);
  }
  END { exit (found_ref && found_nbody && ok) ? 0 : 1; }
' output/zvo_cisajs.dat output/zvo_NBodyG.dat || {
  echo "tJNConserved spin-flip NBodyG does not match cisajs"
  cat output/zvo_cisajs.dat
  cat output/zvo_NBodyG.dat
  exit 1
}

awk -v t="${tol}" '
  NR == FNR {
    if ($1 == 2 && $2 == 1 && $3 == 2 && $4 == 0) {
      ref_re = $5;
      ref_im = $6;
      found_ref = 1;
    }
    next;
  }
  $1 == 1 && $2 == 2 && $3 == 1 && $4 == 2 && $5 == 0 {
    dr = $6 - ref_re;
    di = $7 - ref_im;
    if (dr < 0) dr = -dr;
    if (di < 0) di = -di;
    found_nbody = 1;
    ok = (found_ref && dr < t && di < t);
  }
  END { exit (found_ref && found_nbody && ok) ? 0 : 1; }
' output/zvo_cisajs.dat output/zvo_NBodyG.dat || {
  echo "tJNConserved conjugate spin-flip NBodyG does not match cisajs"
  cat output/zvo_cisajs.dat
  cat output/zvo_NBodyG.dat
  exit 1
}

echo "tJNConserved spin-flip NBodyG matches cisajs."
