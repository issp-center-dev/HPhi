find_program(SPHINX_EXECUTABLE NAMES sphinx-build
  HINTS $ENV{SPHINX_DIR}
  PATH_SUFFIXES bin
  )

if(NOT SPHINX_EXECUTABLE)
  message(FATAL_ERROR "SPHINX_EXECUTABLE (sphinx-build) is not found!")
endif()

set(SPHINX_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/_build)
set(SPHINX_CACHE_DIR ${CMAKE_CURRENT_BINARY_DIR}/_doctrees)
set(SPHINX_HTML_DIR ${CMAKE_CURRENT_BINARY_DIR}/html)
set(SPHINX_PDF_DIR ${CMAKE_CURRENT_BINARY_DIR}/pdf)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/conf.py ${CMAKE_CURRENT_BINARY_DIR})

add_custom_target(doc-en-html ALL
  COMMAND
    ${SPHINX_EXECUTABLE}
    -b html
    -d ${SPHINX_CACHE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${SPHINX_HTML_DIR}
  )

add_custom_target(doc-en-pdf ALL
  COMMAND
    ${SPHINX_EXECUTABLE}
    -b latex
    -d ${SPHINX_CACHE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${SPHINX_PDF_DIR}
  COMMAND
    cd pdf && make
  )

add_custom_target(doc-en DEPENDS doc-en-html doc-en-pdf)

install(DIRECTORY ${SPHINX_HTML_DIR}
  DESTINATION doc
  PATTERN "_sources" EXCLUDE
  PATTERN ".buildinfo" EXCLUDE
  PATTERN "objects.inv" EXCLUDE
  )
install(DIRECTORY ${SPHINX_PDF_DIR}
  DESTINATION doc
  FILES_MATCHING PATTERN "*.pdf"
  )
