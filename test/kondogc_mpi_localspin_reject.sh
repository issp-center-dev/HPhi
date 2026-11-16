#!/bin/sh
set -eu

# KondoGC with a localized spin on an inter-process site must be rejected like
# Kondo; with localized spins on inner sites, CG must match the serial energy.
if [ -z "${MPIRUN}" ]; then
    MPIRUN=""
fi

mkdir -p kondogc_mpi_localspin_reject
cd kondogc_mpi_localspin_reject
cat > stan.in <<INPUT
model = "KondoGC"
method = "CG"
lattice = "chain"
L = 2
t = 1.0
J = 1.0
Lanczos_max = 200
exct = 1
outputmode = "None"
INPUT
../../src/HPhi -sdry stan.in > dry.log 2>&1

write_locspn() {
    {
        echo "================================"
        echo "NlocalSpin 2"
        echo "================================"
        echo "========i_0LocSpn_1IteElc ======"
        echo "================================"
        i=0
        for value in $1; do
            echo "    $i  $value"
            i=$((i + 1))
        done
    } > locspn.def
}

write_locspn "1 1 0 0"
../../src/HPhi -e namelist.def > serial.log 2>&1
cp output/zvo_energy.dat serial_energy.dat

if [ -n "${MPIRUN}" ]; then
    ${MPIRUN} ../../src/HPhi -e namelist.def > mpi.log 2>&1
    paste output/zvo_energy.dat serial_energy.dat |
        awk '$1 == "Energy" {difference=$2-$4; if (difference<0) difference=-difference; if (difference>1e-8) exit 1}'

    write_locspn "0 0 1 1"
    set +e
    ${MPIRUN} ../../src/HPhi -e namelist.def > reject.log 2>&1
    status=$?
    set -e
    test "${status}" -ne 0
    grep -a "local spin in the inter process region" reject.log
fi
