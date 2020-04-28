function(git_clone uri name branch)
    if (NOT IS_DIRECTORY ${name})
        message(STATUS "git clone ${uri}")
        execute_process(COMMAND ${GIT_EXECUTABLE} clone ${uri} ${name}
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
    endif ()
    message(STATUS "git checkout ${branch}")
    execute_process(COMMAND ${GIT_EXECUTABLE} checkout ${branch}
                    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${name}
                   )
endfunction()

git_clone(https://github.com/issp-center-dev/StdFace clone 
         master
         )
