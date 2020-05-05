if(GL_DOWNLOAD_METHOD MATCHES "find|fetch")
    #  Make sure we use GL's own find modules
    list(APPEND CMAKE_MODULE_PATH  ${PROJECT_SOURCE_DIR}/cmake)
    # Make sure find_package looks for the packages we just built instead of rebuilding them
    # when we rerun the cmake config.
    list(APPEND CMAKE_PREFIX_PATH ${CMAKE_INSTALL_PREFIX})

    include(cmake/Fetch_h5pp.cmake)                         # h5pp for writing to file binary in format


    ##################################################################
    ### Link all the things!                                       ###
    ##################################################################
    if(TARGET h5pp::h5pp)
        list(APPEND NATIVE_TARGETS h5pp::h5pp)
    endif()
    if(TARGET openmp::openmp)
        list(APPEND NATIVE_TARGETS openmp::openmp)
    else()
        target_compile_options(project-settings INTERFACE -Wno-unknown-pragmas)
    endif()
    if(TARGET Threads::Threads)
        list(APPEND NATIVE_TARGETS Threads::Threads)
    endif()
    if(NATIVE_TARGETS)
        mark_as_advanced(NATIVE_TARGETS)
    endif()
endif()