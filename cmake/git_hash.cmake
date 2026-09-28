# Write the abbreviated hash of the current commit into a header.
#
# Usage: cmake -DSOURCE_DIR=<top of the source> -DOUTPUT=<header> -P git_hash.cmake
#
# This is run when cmake configures HPhi and every time HPhi is built.
# The header is rewritten only when its content changes, so that nothing is
# recompiled as long as the commit is the same.
# The hash is left empty if the source is not the top of a git repository,
# e.g. a release archive (also one that is extracted inside another repository).

set(HPHI_GIT_HASH "")

find_package(Git QUIET)
if(GIT_FOUND)
  execute_process(
    COMMAND ${GIT_EXECUTABLE} rev-parse --show-toplevel
    WORKING_DIRECTORY ${SOURCE_DIR}
    OUTPUT_VARIABLE _toplevel
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
    RESULT_VARIABLE _rc)
  if(_rc EQUAL 0)
    get_filename_component(_toplevel "${_toplevel}" REALPATH)
    get_filename_component(_source "${SOURCE_DIR}" REALPATH)
    if(_toplevel STREQUAL _source)
      execute_process(
        COMMAND ${GIT_EXECUTABLE} rev-parse --short=8 HEAD
        WORKING_DIRECTORY ${SOURCE_DIR}
        OUTPUT_VARIABLE _hash
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
        RESULT_VARIABLE _rc)
      if(_rc EQUAL 0)
        set(HPHI_GIT_HASH "${_hash}")
      endif()
    endif()
  endif()
endif()

set(_content "#define HPHI_GIT_HASH \"${HPHI_GIT_HASH}\"\n")
set(_old "")
if(EXISTS "${OUTPUT}")
  file(READ "${OUTPUT}" _old)
endif()
if(NOT _content STREQUAL _old)
  file(WRITE "${OUTPUT}" "${_content}")
endif()
