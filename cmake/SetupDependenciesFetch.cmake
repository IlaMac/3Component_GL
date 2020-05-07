if(GL_DOWNLOAD_METHOD MATCHES "find|fetch")
    #  Make sure we use GL's own find modules
    list(APPEND CMAKE_MODULE_PATH  ${PROJECT_SOURCE_DIR}/cmake)
    # Make sure find_package looks for the packages we just built instead of rebuilding them
    # when we rerun the cmake config.
    list(APPEND CMAKE_PREFIX_PATH ${CMAKE_INSTALL_PREFIX})

    include(cmake/Fetch_h5pp.cmake)                         # h5pp for writing to file binary in format


    if(GL_ENABLE_MPI AND NOT TARGET MPI::MPI_CXX)
        find_package(MPI REQUIRED)
    endif()


    ##############################################################################
    ###  OpenMP                                                                ###
    ###  Note that Clang has some  trouble with static openmp and that         ###
    ###  and that static openmp is not recommended. This tries to enable       ###
    ###  static openmp anyway because I find it useful. Installing             ###
    ###  libiomp5 might help for shared linking.                               ###
    ##############################################################################
    if(GL_ENABLE_OPENMP AND NOT TARGET openmp::openmp)
        find_package(OpenMP) # Uses GL's own find module
    endif()


    if(TARGET h5pp::h5pp)
        list(APPEND NATIVE_TARGETS h5pp::h5pp)
    endif()
    if(TARGET Threads::Threads)
        list(APPEND NATIVE_TARGETS Threads::Threads)
    endif()
    if(NATIVE_TARGETS)
        mark_as_advanced(NATIVE_TARGETS)
    endif()
endif()