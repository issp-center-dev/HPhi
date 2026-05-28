#!/bin/sh -e

mkdir -p spinless_onebody_sigma_validation/
cd spinless_onebody_sigma_validation

python3 "$1/test/testSpinlessCalc.py" -p "../../src/HPhi" -m "SpinlessFermion" -s 4

cat > greenone.def <<'EODEF'
===============================
NCisAjs         1
===============================
======== Green functions ======
===============================
0    1    0    0
EODEF

set +e
../../src/HPhi -e namelist.def > run.log 2>&1
rc=$?
set -e

if [ "$rc" -eq 0 ]; then
    echo "FAILED: HPhi succeeded unexpectedly with invalid spin index in OneBodyG"
    cat run.log
    exit 1
fi

if grep -q "OneBodyG spin index must be 0" run.log; then
    echo "PASSED: invalid OneBodyG spin index is rejected for spinless"
    exit 0
fi

echo "FAILED: expected validation error message was not found"
cat run.log
exit 1
