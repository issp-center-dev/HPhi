#!/bin/sh -e

# Test for cmake/git_hash.cmake, which finds the commit printed by "HPhi -v".
#
#   archive (no .git)
#     - cmake/git_archive.txt filled in with a hash -> its first 8 digits
#     - cmake/git_archive.txt left as "$Format:%H$"  -> empty
#     - no cmake/git_archive.txt                     -> empty
#   git repository
#     - clean                                        -> first 8 digits of HEAD
#     - a file under the version control is changed  -> followed by "-dirty"
#     - a file which is not under the version control is added -> not dirty
#     - the header is left untouched when its content does not change
#
# Usage: version_git_hash.sh <top of the source> <cmake>

testname="version_git_hash"
srcdir=$1
cmake=${2:-cmake}
script="${srcdir}/cmake/git_hash.cmake"

rm -rf "${testname}"
mkdir -p "${testname}"
cd "${testname}"
workdir="$(pwd)"

fail() {
  echo "FAILED (${testname}): $1" >&2
  exit 1
}

# expect <label> <source directory> <expected value of HPHI_GIT_HASH>
expect() {
  "${cmake}" -DSOURCE_DIR="$2" -DOUTPUT="${workdir}/version_git.h" -P "${script}" \
    || fail "$1: cmake failed"
  actual=$(cat "${workdir}/version_git.h")
  [ "x${actual}" = "x#define HPHI_GIT_HASH \"$3\"" ] || \
    fail "$1: '${actual}', expected HPHI_GIT_HASH \"$3\""
  echo "$1: ${actual}"
}

#
# Archives
#
mkdir -p archive/cmake
echo "0123456789abcdef0123456789abcdef01234567" > archive/cmake/git_archive.txt
expect "archive, filled in" "${workdir}/archive" "01234567"

echo '$Format:%H$' > archive/cmake/git_archive.txt
expect "archive, not filled in" "${workdir}/archive" ""

rm archive/cmake/git_archive.txt
expect "no git_archive.txt" "${workdir}/archive" ""

# The file in the source has to be the one which git archive fills in
[ "x$(cat "${srcdir}/cmake/git_archive.txt")" = 'x$Format:%H$' ] || {
  # in an archive it is filled in with a hash
  grep -q -E '^[0-9a-f]{40}$' "${srcdir}/cmake/git_archive.txt" || \
    fail "cmake/git_archive.txt is neither \$Format:%H\$ nor a hash"
}
if [ -f "${srcdir}/.gitattributes" ]; then
  grep -q -E '^cmake/git_archive.txt[[:space:]]+export-subst' "${srcdir}/.gitattributes" || \
    fail "export-subst is not set to cmake/git_archive.txt in .gitattributes"
fi

#
# Git repository
#
if ! git --version > /dev/null 2>&1; then
  echo "git is not available. The cases of a git repository are skipped."
  exit 0
fi

mkdir -p repo
(
  cd repo
  git init -q .
  echo "a" > tracked.txt
  git add tracked.txt
  git -c user.name=test -c user.email=test@example.com -c commit.gpgsign=false \
    commit -q -m "test"
) || fail "a git repository could not be made"
hash=$(cd repo && git rev-parse HEAD | cut -c 1-8)

expect "repository, clean" "${workdir}/repo" "${hash}"

# The header is not rewritten when the content is the same
cp -p version_git.h version_git.h.keep
sleep 1
expect "repository, clean (again)" "${workdir}/repo" "${hash}"
[ version_git.h -nt version_git.h.keep ] && fail "the header was rewritten without a change"

echo "b" > repo/untracked.txt
expect "repository, untracked file" "${workdir}/repo" "${hash}"

echo "c" >> repo/tracked.txt
expect "repository, changed" "${workdir}/repo" "${hash}-dirty"

echo "The commit is found as expected."
