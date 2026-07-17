#!/bin/sh
set -eu

testname="modpara_parse_validation"
hphi="../../src/HPhi"

rm -rf "${testname}"
mkdir -p "${testname}"
cd "${testname}"

cat > stan.in <<EOF
model = "Hubbard"
method = "Lanczos"
lattice = "chain"
L = 4
t = 1.0
U = 4.0
nelec = 4
2Sz = 0
EOF

"${hphi}" -sdry stan.in > generate.log 2>&1 || {
  cat generate.log
  exit 1
}

awk '
  $1 == "Lanczos_max" {
    print "Lanczos_max"
    next
  }
  { print }
' modpara.def > modpara.def.tmp
mv modpara.def.tmp modpara.def

set +e
"${hphi}" -e namelist.def > malformed.log 2>&1
rc=$?
set -e

if [ "${rc}" -eq 0 ]; then
  echo "Expected a ModPara line without a numeric value to be rejected."
  cat malformed.log
  exit 1
fi

if ! grep -q "ModPara line must contain a keyword and numeric value" malformed.log; then
  echo "Malformed ModPara was not rejected by the expected parser guard."
  cat malformed.log
  exit 1
fi

echo "ModPara parser rejects a line without a numeric value."
