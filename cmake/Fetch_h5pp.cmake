if(NOT TARGET h5pp::h5pp AND GL_DOWNLOAD_METHOD STREQUAL "find")
    find_package(h5pp 1.7.1 QUIET
            HINTS ${CMAKE_INSTALL_PREFIX}
            PATH_SUFFIXES h5pp
            )
    if(h5pp_FOUND AND TARGET h5pp::h5pp)
        message(STATUS "Found h5pp")
    endif()
endif()

if(NOT TARGET h5pp::h5pp AND GL_DOWNLOAD_METHOD MATCHES "fetch")
    find_package(h5pp 1.7.1
            HINTS ${CMAKE_INSTALL_PREFIX}
            PATH_SUFFIXES h5pp
            NO_CMAKE_PACKAGE_REGISTRY)
    if(h5pp_FOUND AND TARGET h5pp::h5pp)
        message(STATUS "Found h5pp")
    endif()
endif()

if(NOT TARGET h5pp::h5pp AND GL_DOWNLOAD_METHOD MATCHES "fetch")
    message(STATUS "h5pp will be installed into ${CMAKE_INSTALL_PREFIX}")
    include(${PROJECT_SOURCE_DIR}/cmake/BuildDependency.cmake)
    list(APPEND H5PP_CMAKE_OPTIONS  -DCMAKE_PREFIX_PATH:PATH=${CMAKE_PREFIX_PATH})
    list(APPEND H5PP_CMAKE_OPTIONS  -DH5PP_ENABLE_EIGEN3:PATH=${GL_ENABLE_EIGEN3})
    list(APPEND H5PP_CMAKE_OPTIONS  -DH5PP_ENABLE_SPDLOG:PATH=${GL_ENABLE_SPDLOG})
    list(APPEND H5PP_CMAKE_OPTIONS  -DH5PP_DOWNLOAD_METHOD:BOOL=${GL_DOWNLOAD_METHOD})
    list(APPEND H5PP_CMAKE_OPTIONS  -DH5PP_PRINT_INFO:PATH=${GL_PRINT_INFO})
    list(APPEND H5PP_CMAKE_OPTIONS  -DH5PP_PREFER_CONDA_LIBS:BOOL=${GL_PREFER_CONDA_LIBS})
    build_dependency(h5pp "${CMAKE_INSTALL_PREFIX}" "${H5PP_CMAKE_OPTIONS}")

    find_package(h5pp 1.7.1 HINTS
            HINTS ${CMAKE_INSTALL_PREFIX}
            PATH_SUFFIXES h5pp
            NO_DEFAULT_PATH)
    if(h5pp_FOUND AND TARGET h5pp::h5pp)
        message(STATUS "h5pp installed successfully")
    else()
        message(FATAL_ERROR "h5pp could not be installed")
    endif()
endif()
