#!/bin/sh -e

# Test for the last part of dist.sh, which makes the release archive and
# writes the hash of the commit into cmake/git_archive.txt of the archive.
#
#   - A failure of the archive creation stops dist.sh with the status of
#     git-archive-all, and the following steps are not run.
#   - On success the hash of HEAD is written into the archive.
#   - A failure of the repacking stops dist.sh with the status of tar, and the
#     archive made by git-archive-all is left as it was.
#   - An archive which is filled in already is not touched.
#   - The temporary files are removed in every case.
#
# dist.sh is run in a small tree made for the test, with "make" and
# "git-archive-all" replaced by stubs, so that neither the manuals nor
# git-archive-all are needed.
#
# Usage: dist_archive_failure.sh <top of the source>

testname="dist_archive_failure"
srcdir=$1

rm -rf "${testname}"
mkdir -p "${testname}"
cd "${testname}"
workdir="$(pwd)"

fail() {
  echo "FAILED (${testname}): $1" >&2
  exit 1
}

if ! git --version > /dev/null 2>&1; then
  echo "git is not available. The test is skipped."
  exit 0
fi

#
# A small tree with what dist.sh reads
#
mkdir -p root/doc/ja root/doc/en root/doc/tutorial/en root/src/include root/cmake
cp "${srcdir}/dist.sh" root/
cp "${srcdir}/src/include/version.h" root/src/include/
echo '$Format:%H$' > root/cmake/git_archive.txt
(
  cd root
  git init -q .
  git add dist.sh src/include/version.h cmake/git_archive.txt
  git -c user.name=test -c user.email=test@example.com -c commit.gpgsign=false \
    commit -q -m "test"
) || fail "a git repository could not be made"
hash=$(cd root && git rev-parse HEAD)
major=$(awk '$2=="HPHI_VERSION_MAJOR"{print $3}' root/src/include/version.h)
minor=$(awk '$2=="HPHI_VERSION_MINOR"{print $3}' root/src/include/version.h)
patch=$(awk '$2=="HPHI_VERSION_PATCH"{print $3}' root/src/include/version.h)
tarball="HPhi-${major}.${minor}.${patch}.tar.gz"
prefix="HPhi-${major}.${minor}.${patch}"

#
# Stubs. STUB_ARCHIVE and STUB_TAR choose how they behave.
#
mkdir -p stub
realtar=$(command -v tar)
cat > stub/make <<STUB
#!/bin/sh
exit 0
STUB
cat > stub/git-archive-all <<STUB
#!/bin/sh
[ "\$1" = "--help" ] && exit 0
[ "\${STUB_ARCHIVE}" = "fail" ] && exit 23
for last in "\$@"; do :; done
rm -rf "${workdir}/stage"
mkdir -p "${workdir}/stage/${prefix}/cmake"
if [ "\${STUB_ARCHIVE}" = "filled" ]; then
  echo "0123456789abcdef0123456789abcdef01234567" > "${workdir}/stage/${prefix}/cmake/git_archive.txt"
elif [ "\${STUB_ARCHIVE}" != "nofile" ]; then
  cp cmake/git_archive.txt "${workdir}/stage/${prefix}/cmake/"
fi
"${realtar}" czf "\${last}" -C "${workdir}/stage" "${prefix}"
STUB
cat > stub/tar <<STUB
#!/bin/sh
if [ "\${STUB_TAR}" = "fail_create" ] && [ "\$1" = "czf" ]; then
  exit 5
fi
exec "${realtar}" "\$@"
STUB
chmod +x stub/make stub/git-archive-all stub/tar

# run_dist <STUB_ARCHIVE> <STUB_TAR>; the exit status is left in ${rc}
run_dist() {
  rm -f "root/${tarball}" "root/${tarball}.tmp"
  set +e
  (
    cd root
    TMPDIR="${workdir}/tmp" STUB_ARCHIVE="$1" STUB_TAR="$2" PATH="${workdir}/stub:${PATH}" \
      sh dist.sh > "${workdir}/dist.log" 2>&1
  )
  rc=$?
  set -e
}

# content of cmake/git_archive.txt in the archive
archive_txt() {
  "${realtar}" xzf "root/${tarball}" -O "${prefix}/cmake/git_archive.txt"
}

# no temporary file is left
check_clean() {
  [ ! -e "root/${tarball}.tmp" ] || fail "$1: ${tarball}.tmp is left"
  [ -z "$(ls -A "${workdir}/tmp")" ] || fail "$1: a temporary directory is left"
}

mkdir -p tmp

label="archive creation fails"
run_dist fail none
[ "${rc}" = "23" ] || { cat dist.log; fail "${label}: exit status ${rc}, expected 23"; }
[ ! -e "root/${tarball}" ] || fail "${label}: an archive exists"
check_clean "${label}"
echo "${label}: exit status ${rc}"

label="success"
run_dist ok none
[ "${rc}" = "0" ] || { cat dist.log; fail "${label}: exit status ${rc}"; }
[ "x$(archive_txt)" = "x${hash}" ] || fail "${label}: git_archive.txt is '$(archive_txt)', expected ${hash}"
check_clean "${label}"
echo "${label}: the hash is written"

label="repacking fails"
run_dist ok fail_create
[ "${rc}" = "5" ] || { cat dist.log; fail "${label}: exit status ${rc}, expected 5"; }
[ "x$(archive_txt)" = 'x$Format:%H$' ] || fail "${label}: the archive was changed"
check_clean "${label}"
echo "${label}: exit status ${rc}, the archive is left as it was"

label="filled in already"
run_dist filled fail_create
[ "${rc}" = "0" ] || { cat dist.log; fail "${label}: exit status ${rc}"; }
[ "x$(archive_txt)" = "x0123456789abcdef0123456789abcdef01234567" ] || fail "${label}: the archive was changed"
check_clean "${label}"
echo "${label}: not touched"

label="no git_archive.txt in the archive"
run_dist nofile none
[ "${rc}" != "0" ] || { cat dist.log; fail "${label}: exit status 0"; }
check_clean "${label}"
echo "${label}: exit status ${rc}"

echo "dist.sh propagates the failures."
