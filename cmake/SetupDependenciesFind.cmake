if(GL_DOWNLOAD_METHOD MATCHES "find")
    # Let cmake find our Find<package>.cmake modules
    list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
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



if(GL_ENABLE_MPI)
    if(NOT BUILD_SHARED_LIBS)
        message(WARNING "Linking MPI statically is discouraged and may fail. Try --enable-shared or setting BUILD_SHARED_LIBS=ON.")
    endif()
    find_package(MPI REQUIRED)
    list(APPEND NATIVE_TARGETS MPI::MPI_CXX)
endif()


##############################################################################
###  OpenMP                                                                ###
###  Note that Clang has some  trouble with static openmp and that         ###
###  and that static openmp is not recommended. This tries to enable       ###
###  static openmp anyway because I find it useful. Installing             ###
###  libiomp5 might help for shared linking.                               ###
##############################################################################
if(GL_ENABLE_OPENMP)
    find_package(OpenMP) # Uses GL's own find module
    list(APPEND NATIVE_TARGETS openmp::openmp)
endif()