if(GL_DOWNLOAD_METHOD MATCHES "find")
    # Let cmake find our Find<package>.cmake modules
    list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
    list(APPEND CMAKE_PREFIX_PATH ${CMAKE_INSTALL_PREFIX})

    if(GL_PREFER_CONDA_LIBS)
        list(APPEND CMAKE_PREFIX_PATH
                $ENV{CONDA_PREFIX}
                $ENV{HOME}/anaconda3
                $ENV{HOME}/anaconda
                $ENV{HOME}/miniconda3
                $ENV{HOME}/miniconda
                $ENV{HOME}/.conda
                )
        endif()

    list(APPEND CMAKE_PREFIX_PATH
            $ENV{EBROOTOPENMPI}
            $ENV{EBROOTHDF5}
            $ENV{EBROOTEIGEN}
            $ENV{EBROOTSPDLOG}
            )
endif()


# Things below  should happen regardless of download mode find|fetch|find-or-fetch|conan


if(GL_ENABLE_MPI AND NOT TARGET MPI::MPI_CXX)
    find_package(MPI REQUIRED)
    if(TARGET MPI::MPI_CXX)
        list(APPEND NATIVE_TARGETS MPI::MPI_CXX)
    endif()
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
    if(TARGET openmp::openmp)
        list(APPEND NATIVE_TARGETS openmp::openmp)
    else()
        target_compile_options(project-settings INTERFACE -Wno-unknown-pragmas)
    endif()
endif()