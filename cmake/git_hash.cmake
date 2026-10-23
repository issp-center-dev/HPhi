# Write the commit which HPhi is built from into a header.
#
# Usage: cmake -DSOURCE_DIR=<top of the source> -DOUTPUT=<header> -P git_hash.cmake
#
# The header defines HPHI_GIT_HASH as the first 8 digits of the hash of
#   - HEAD, in a git repository
#   - the commit which the archive was made from, in an archive
#     (cmake/git_archive.txt is filled in by "git archive", see .gitattributes,
#     and by dist.sh)
# and as an empty string otherwise.
# In a git repository "-dirty" follows the hash if files under the version
# control have changes which are not committed, as "211bd466-dirty". Files
# which are not under the version control are not taken into account, as
# "git describe --dirty".
#
# This is run when cmake configures HPhi and every time HPhi is built.
# The header is rewritten only when its content changes, so that nothing is
# recompiled as long as the commit is the same.

set(_hash "")
set(_dirty "")

if(EXISTS "${SOURCE_DIR}/.git")
  find_package(Git QUIET)
  if(GIT_FOUND)
    execute_process(
      COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
      WORKING_DIRECTORY "${SOURCE_DIR}"
      OUTPUT_VARIABLE _hash
      OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_QUIET
      RESULT_VARIABLE _rc)
    if(NOT _rc EQUAL 0)
      set(_hash "")
    endif()
    execute_process(
      COMMAND ${GIT_EXECUTABLE} status --porcelain --untracked-files=no
              --ignore-submodules=untracked
      WORKING_DIRECTORY "${SOURCE_DIR}"
      OUTPUT_VARIABLE _status
      OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_QUIET
      RESULT_VARIABLE _rc)
    if(_rc EQUAL 0 AND NOT "${_status}" STREQUAL "")
      set(_dirty "-dirty")
    endif()
  endif()
elseif(EXISTS "${SOURCE_DIR}/cmake/git_archive.txt")
  file(STRINGS "${SOURCE_DIR}/cmake/git_archive.txt" _hash LIMIT_COUNT 1)
  string(STRIP "${_hash}" _hash)
endif()

# "$Format:%H$" is left in git_archive.txt if it is not filled in
string(LENGTH "${_hash}" _length)
if(_hash MATCHES "^[0-9a-f]+$" AND NOT _length LESS 8)
  string(SUBSTRING "${_hash}" 0 8 _hash)
  set(HPHI_GIT_HASH "${_hash}${_dirty}")
else()
  set(HPHI_GIT_HASH "")
endif()

set(_content "#define HPHI_GIT_HASH \"${HPHI_GIT_HASH}\"\n")
set(_old "")
if(EXISTS "${OUTPUT}")
  file(READ "${OUTPUT}" _old)
endif()
if(NOT _content STREQUAL _old)
  file(WRITE "${OUTPUT}" "${_content}")
endif()
