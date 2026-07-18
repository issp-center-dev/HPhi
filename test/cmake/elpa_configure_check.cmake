if(NOT DEFINED HPHI_SOURCE_DIR OR
   NOT DEFINED HPHI_TEST_ROOT OR
   NOT DEFINED HPHI_PRELOAD_CACHE)
  message(FATAL_ERROR "Missing ELPA configure-test arguments")
endif()

get_filename_component(_test_root_name "${HPHI_TEST_ROOT}" NAME)
if(NOT _test_root_name STREQUAL "elpa_configure_check")
  message(FATAL_ERROR
    "Refusing to remove an unexpected ELPA configure-test root: ${HPHI_TEST_ROOT}")
endif()

file(REMOVE_RECURSE "${HPHI_TEST_ROOT}")
file(MAKE_DIRECTORY "${HPHI_TEST_ROOT}")
include("${HPHI_PRELOAD_CACHE}")

# The parent build may expose more than one ELPA include directory. The
# negative-library probe accepts a single ELPA_INCLUDE_DIR cache value, so
# select the actual API directory rather than forwarding a semicolon list as
# multiple execute_process arguments.
set(_valid_elpa_include_dir)
foreach(_elpa_include_dir ${HPHI_VALID_ELPA_INCLUDE_DIRS})
  if(EXISTS "${_elpa_include_dir}/elpa/elpa.h")
    set(_valid_elpa_include_dir "${_elpa_include_dir}")
    break()
  endif()
endforeach()
if(NOT _valid_elpa_include_dir)
  message(FATAL_ERROR
    "No ELPA include directory containing elpa/elpa.h was provided by the parent build")
endif()

function(hphi_expect_configure_failure name expected_text)
  set(_build "${HPHI_TEST_ROOT}/${name}")
  file(MAKE_DIRECTORY "${_build}")
  execute_process(
    COMMAND "${CMAKE_COMMAND}" -C "${HPHI_PRELOAD_CACHE}"
            -DGIT_SUBMODULE_UPDATE=OFF ${ARGN} "${HPHI_SOURCE_DIR}"
    WORKING_DIRECTORY "${_build}"
    RESULT_VARIABLE _result
    OUTPUT_VARIABLE _stdout
    ERROR_VARIABLE _stderr)
  set(_log "${_stdout}\n${_stderr}")
  if(_result EQUAL 0)
    message(FATAL_ERROR
      "${name}: configure unexpectedly succeeded\n${_log}")
  endif()
  string(FIND "${_log}" "${expected_text}" _message_pos)
  if(_message_pos EQUAL -1)
    message(FATAL_ERROR
      "${name}: expected diagnostic not found: ${expected_text}\n${_log}")
  endif()
endfunction()

hphi_expect_configure_failure(
  explicit_scalapack_off
  "USE_ELPA=ON requires ScaLAPACK, but USE_SCALAPACK=OFF was set explicitly"
  -DUSE_ELPA=ON -DUSE_SCALAPACK=OFF)

hphi_expect_configure_failure(
  config_elpa_with_explicit_scalapack_off
  "USE_ELPA=ON requires ScaLAPACK, but USE_SCALAPACK=OFF was set explicitly"
  -DCONFIG=elpa -DUSE_SCALAPACK=OFF)

hphi_expect_configure_failure(
  missing_mpi
  "USE_ELPA=ON requires a working MPI C implementation"
  -DUSE_ELPA=ON -DCMAKE_DISABLE_FIND_PACKAGE_MPI=TRUE)

hphi_expect_configure_failure(
  invalid_scalapack
  "Could NOT find ScaLAPACK"
  -DUSE_ELPA=ON -DSCALAPACK_LIBRARIES=-lm)

hphi_expect_configure_failure(
  invalid_elpa_paths
  "ELPA_INCLUDE_DIR does not contain elpa/elpa.h"
  -DUSE_ELPA=ON
  -DELPA_INCLUDE_DIR=/hphi/nonexistent/elpa/include
  -DELPA_LIBRARY=/hphi/nonexistent/lib/libelpa.so)

# An existing file is not sufficient: the base elpa_init API must link.
set(_fake_elpa_library "${HPHI_TEST_ROOT}/libelpa-invalid.so")
file(WRITE "${_fake_elpa_library}" "not an ELPA library\n")
hphi_expect_configure_failure(
  invalid_elpa_library
  "ELPA was found, but a test program using elpa_init could not be linked"
  -DUSE_ELPA=ON
  "-DELPA_INCLUDE_DIR=${_valid_elpa_include_dir}"
  "-DELPA_LIBRARY=${_fake_elpa_library}")

# A prior explicit OFF is sticky so a later bare -DUSE_ELPA=ON reports a
# conflict. A configuration file that deliberately enables ScaLAPACK must,
# however, clear that stale marker and configure successfully in the same
# build directory.
set(_sticky_build "${HPHI_TEST_ROOT}/config_reenable")
file(MAKE_DIRECTORY "${_sticky_build}")
execute_process(
  COMMAND "${CMAKE_COMMAND}" -C "${HPHI_PRELOAD_CACHE}"
          -DGIT_SUBMODULE_UPDATE=OFF -DUSE_SCALAPACK=OFF "${HPHI_SOURCE_DIR}"
  WORKING_DIRECTORY "${_sticky_build}"
  RESULT_VARIABLE _first_result
  OUTPUT_VARIABLE _first_stdout
  ERROR_VARIABLE _first_stderr)
if(NOT _first_result EQUAL 0)
  message(FATAL_ERROR
    "config_reenable: initial configure failed\n${_first_stdout}\n${_first_stderr}")
endif()

execute_process(
  COMMAND "${CMAKE_COMMAND}" -C "${HPHI_PRELOAD_CACHE}"
          -DGIT_SUBMODULE_UPDATE=OFF -DCONFIG=elpa "${HPHI_SOURCE_DIR}"
  WORKING_DIRECTORY "${_sticky_build}"
  RESULT_VARIABLE _second_result
  OUTPUT_VARIABLE _second_stdout
  ERROR_VARIABLE _second_stderr)
if(NOT _second_result EQUAL 0)
  message(FATAL_ERROR
    "config_reenable: CONFIG=elpa did not clear the prior explicit-OFF marker\n"
    "${_second_stdout}\n${_second_stderr}")
endif()

file(REMOVE_RECURSE "${HPHI_TEST_ROOT}")
