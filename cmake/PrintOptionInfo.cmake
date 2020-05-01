if(GL_PRINT_INFO)
    # Print CMake options
    message(STATUS  "|----------------\n"
            "-- | BUILD_SHARED_LIBS       : ${BUILD_SHARED_LIBS}\n"
            "-- | CMAKE_BUILD_TYPE        : ${CMAKE_BUILD_TYPE}\n"
            "-- | CMAKE_INSTALL_PREFIX    : ${CMAKE_INSTALL_PREFIX}\n"
            "-- | CMAKE_PREFIX_PATH       : ${CMAKE_PREFIX_PATH}\n"
            "-- | GL_ENABLE_TESTS         : ${GL_ENABLE_TESTS}\n"
            "-- | GL_ENABLE_H5PP          : ${GL_ENABLE_H5PP}\n"
            "-- | GL_ENABLE_EIGEN3        : ${GL_ENABLE_EIGEN3}\n"
            "-- | GL_ENABLE_SPDLOG        : ${GL_ENABLE_SPDLOG}\n"
            "-- | GL_ENABLE_MPI           : ${GL_ENABLE_MPI}\n"
            "-- | GL_ENABLE_OPENMP        : ${GL_ENABLE_OPENMP}\n"
            "-- | GL_ENABLE_LTO           : ${GL_ENABLE_LTO}\n"
            "-- | GL_DOWNLOAD_METHOD      : ${GL_DOWNLOAD_METHOD}\n"
            "-- | GL_PREFER_CONDA_LIBS    : ${GL_PREFER_CONDA_LIBS}\n"
            "-- | GL_PRINT_INFO           : ${GL_PRINT_INFO}\n")
endif ()
