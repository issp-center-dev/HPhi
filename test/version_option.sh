#!/bin/sh -e

# Test for "HPhi -v" and "HPhi --version".
#
#   - The exit status is 0, so that the option can be used as a probe.
#   - The output is "HPhi version X.Y.Z" or "HPhi version X.Y.Z (hash)",
#     where X.Y.Z agrees with src/include/version.h and the hash is the
#     abbreviated hash (8 characters or more) of a git commit.
#
# The hash is not compared with that of the source tree: it is the one at the
# time of the build, and the tests may be run after another commit.
#
# Usage: version_option.sh <top of the source>

testname="version_option"
srcdir=$1
header="${srcdir}/src/include/version.h"

mkdir -p "${testname}"
cd "${testname}"

fail() {
  echo "FAILED (${testname}): $1" >&2
  exit 1
}

[ -f "${header}" ] || fail "${header} is not found"
major=$(awk '$2=="HPHI_VERSION_MAJOR"{print $3}' "${header}")
minor=$(awk '$2=="HPHI_VERSION_MINOR"{print $3}' "${header}")
patch=$(awk '$2=="HPHI_VERSION_PATCH"{print $3}' "${header}")
expected="HPhi version ${major}.${minor}.${patch}"

for option in -v --version; do
  set +e
  ../../src/HPhi ${option} > version.log 2>&1
  rc=$?
  set -e
  [ "${rc}" = "0" ] || { cat version.log; fail "${option}: exit status ${rc}"; }

  line=$(grep "HPhi version" version.log | head -1)
  echo "${option}: ${line}"

  number=$(printf '%s\n' "${line}" | sed 's/ *(.*$//; s/ *$//')
  [ "x${number}" = "x${expected}" ] || fail "${option}: '${number}' != '${expected}'"

  case "${line}" in
    *"("*)
      printf '%s\n' "${line}" | grep -q -E '^HPhi version [0-9]+\.[0-9]+\.[0-9]+ \([0-9a-f]{8,}\)$' || \
        fail "${option}: the hash is not in the form (xxxxxxxx)"
      ;;
  esac
done

# Without an option HPhi prints the usage, which is an error
set +e
../../src/HPhi > usage.log 2>&1
rc=$?
set -e
[ "${rc}" != "0" ] || fail "no option: exit status 0"

echo "The version is printed as expected."
