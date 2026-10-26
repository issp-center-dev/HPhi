#!/bin/sh -e

# Regression test for Kondo + Ncond without 2Sz (KondoNConserved) with
# local spins placed at arbitrary sites.
#
# The basis offsets of KondoNConserved used to be right only for the layout
# written by the standard mode (all local spins in the first half of the sites,
# all itinerant sites in the second half).  For any other layout the basis had
# the right dimension but overlapping entries, so the spectrum was wrong while
# kondo_nconserved_dimension still passed.
#
# For each layout the full spectrum of KondoNConserved is compared with the
# union of the spectra of the canonical Kondo model over all 2Sz sectors, which
# is built by an independent basis routine.  The first case is in addition
# compared with a ground-state energy obtained by an independent exact
# diagonalization.  The test is serial-only: it verifies basis construction,
# not MPI decomposition.

testname="fulldiag_kondo_nconserved_layout"
tolerance="0.00000001"
hphi="$(pwd)/../src/HPhi"

mkdir -p "${testname}"
cd "${testname}"
workdir="$(pwd)"

fail() {
  echo "FAILED (${testname}): $1" >&2
  exit 1
}

# write_input <layout> <ncond> <2Sz or "none"> <uniform>
#   layout : string of 0 (itinerant) and 1 (local spin), one character per site
#   uniform: 1 -> t = 2 and J = 1 on all bonds, no coupling between local spins
#            0 -> bond-dependent t and J, and couplings between local spins
write_input() {
  layout=$1
  ncond=$2
  twosz=$3
  uniform=$4

  awk -v layout="${layout}" -v ncond="${ncond}" -v twosz="${twosz}" -v uniform="${uniform}" '
  function abs(x) { return x < 0 ? -x : x }
  function add_heisenberg(a, b, J,    s, sp) {
    for (s = 0; s < 2; s++) {
      for (sp = 0; sp < 2; sp++) {
        inter[nint++] = sprintf("%d %d %d %d %d %d %d %d %.10f 0.0", a, s, a, s, b, sp, b, sp, (s == sp ? J : -J) / 4)
      }
    }
    inter[nint++] = sprintf("%d 0 %d 1 %d 1 %d 0 %.10f 0.0", a, a, b, b, J / 2)
    inter[nint++] = sprintf("%d 0 %d 1 %d 1 %d 0 %.10f 0.0", b, b, a, a, J / 2)
  }
  BEGIN {
    nsite = length(layout)
    nitin = 0
    nloc = 0
    for (i = 0; i < nsite; i++) {
      loc[i] = substr(layout, i + 1, 1)
      if (loc[i] == 0) itin[nitin++] = i
      else locs[nloc++] = i
    }

    f = "locspn.def"
    printf "====\nNlocalSpin %d\n====\n====\n====\n", nloc > f
    for (i = 0; i < nsite; i++) printf "%d %d\n", i, loc[i] > f
    close(f)

    ntrans = 0
    for (k = 0; k < nitin - 1; k++) {
      t = (uniform == 1 ? 2.0 : 1.0 + 0.1 * k)
      for (s = 0; s < 2; s++) {
        trans[ntrans++] = sprintf("%d %d %d %d %.10f 0.0", itin[k], s, itin[k + 1], s, t)
        trans[ntrans++] = sprintf("%d %d %d %d %.10f 0.0", itin[k + 1], s, itin[k], s, t)
      }
    }
    f = "trans.def"
    printf "====\nNTransfer %d\n====\n====\n====\n", ntrans > f
    for (i = 0; i < ntrans; i++) print trans[i] > f
    close(f)

    nint = 0
    for (k = 0; k < nloc; k++) {
      # couple each local spin to the nearest itinerant site
      best = itin[0]
      for (m = 1; m < nitin; m++) {
        if (abs(itin[m] - locs[k]) < abs(best - locs[k])) best = itin[m]
      }
      add_heisenberg(locs[k], best, (uniform == 1 ? 1.0 : 1.0 + 0.17 * k))
    }
    if (uniform != 1) {
      for (k = 0; k < nloc - 1; k++) add_heisenberg(locs[k], locs[k + 1], 0.3 + 0.05 * k)
    }
    f = "interall.def"
    printf "====\nNInterAll %d\n====\n====\n====\n", nint > f
    for (i = 0; i < nint; i++) print inter[i] > f
    close(f)

    f = "modpara.def"
    printf "----\nModel_Parameters 0\n----\nHPhi_Cal_Parameters\n----\n" > f
    printf "CDataFileHead zvo\nCParaFileHead zqp\n----\n" > f
    printf "Nsite %d\nNcond %d\n", nsite, ncond > f
    if (twosz != "none") printf "2Sz %d\n", twosz > f
    printf "Lanczos_max 2000\ninitial_iv -1\nexct 1\nLanczosEps 14\nLanczosTarget 2\n" > f
    close(f)

    f = "calcmod.def"
    printf "CalcType 2\nCalcModel 2\n" > f
    close(f)

    f = "namelist.def"
    printf "ModPara modpara.def\nLocSpin locspn.def\nTrans trans.def\n" > f
    printf "InterAll interall.def\nCalcMod calcmod.def\n" > f
    close(f)
  }'
}

# run_hphi <label> <directory> <layout> <ncond> <2Sz or "none"> <uniform>
run_hphi() {
  run_label=$1
  run_dir=$2

  rm -rf "${run_dir}"
  mkdir -p "${run_dir}"
  (
    cd "${run_dir}"
    write_input "$3" "$4" "$5" "$6"
    "${hphi}" -e namelist.def > run.log 2>&1 || { cat run.log; fail "${run_label}: HPhi failed"; }
    [ -f output/Eigenvalue.dat ] || { cat run.log; fail "${run_label}: Eigenvalue.dat was not written"; }
  )
}

# run_case <layout> <ncond> <uniform> [reference ground-state energy]
run_case() {
  layout=$1
  ncond=$2
  uniform=$3
  reference=$4
  label="layout_${layout}_ncond${ncond}"
  casedir="${workdir}/${label}"

  rm -rf "${casedir}"
  mkdir -p "${casedir}"

  nloc=$(printf '%s' "${layout}" | tr -d '0' | wc -c | tr -d ' ')
  ne=$((ncond + nloc))

  # Reference: canonical Kondo in every 2Sz sector
  : > "${casedir}/sectors.dat"
  twosz=$((-ne))
  while [ "${twosz}" -le "${ne}" ]; do
    run_hphi "${label} 2Sz=${twosz}" "${casedir}/sz_${twosz}" "${layout}" "${ncond}" "${twosz}" "${uniform}"
    awk 'NF >= 2 { print $2 }' "${casedir}/sz_${twosz}/output/Eigenvalue.dat" >> "${casedir}/sectors.dat"
    twosz=$((twosz + 2))
  done
  sort -n "${casedir}/sectors.dat" > "${casedir}/expected.dat"

  # Ncond without 2Sz: promoted to KondoNConserved
  run_hphi "${label} no 2Sz" "${casedir}/nconserved" "${layout}" "${ncond}" none "${uniform}"
  awk 'NF >= 2 { print $2 }' "${casedir}/nconserved/output/Eigenvalue.dat" | sort -n > "${casedir}/actual.dat"

  nexpected=$(wc -l < "${casedir}/expected.dat" | tr -d ' ')
  nactual=$(wc -l < "${casedir}/actual.dat" | tr -d ' ')
  [ "${nexpected}" -gt 0 ] || fail "${label}: the reference spectrum is empty"
  [ "x${nactual}" = "x${nexpected}" ] || \
    fail "${label}: ${nactual} eigenvalues, expected ${nexpected}"

  paste "${casedir}/actual.dat" "${casedir}/expected.dat" > "${casedir}/paste.dat"
  maxdiff=$(awk 'BEGIN { m = 0 } { d = $1 - $2; if (d < 0) d = -d; if (d > m) m = d } END { printf "%.12f", m }' "${casedir}/paste.dat")
  ok=$(awk -v d="${maxdiff}" -v tol="${tolerance}" 'BEGIN { print (d < tol) ? 1 : 0 }')
  [ "${ok}" -eq 1 ] || {
    cat "${casedir}/paste.dat"
    fail "${label}: spectrum differs from the canonical sectors (max diff ${maxdiff})"
  }

  if [ -n "${reference}" ]; then
    e0=$(head -1 "${casedir}/actual.dat")
    ok=$(awk -v a="${e0}" -v b="${reference}" -v tol="${tolerance}" 'BEGIN { d = a - b; if (d < 0) d = -d; print (d < tol) ? 1 : 0 }')
    [ "${ok}" -eq 1 ] || fail "${label}: ground-state energy ${e0} != reference ${reference}"
  fi

  echo "${label}: ${nactual} eigenvalues agree (max diff ${maxdiff})"
}

# Alternating layout, t = 2, J = 1: reference from an independent exact diagonalization
run_case 0101 2 1 -4.1219214771
# Layout of the standard mode (worked before the fix as well)
run_case 1100 2 0
# Local spins in the second half
run_case 0011 2 0
# Odd number of sites: the central site is itinerant / a local spin
run_case 10010 2 0
run_case 00101 3 0
# Dilute local spins, itinerant sites in both halves
run_case 110000 3 0
run_case 100010 3 0

echo "KondoNConserved spectra agree with the canonical Kondo sectors."
