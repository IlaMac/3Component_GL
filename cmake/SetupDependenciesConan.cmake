
if(GL_DOWNLOAD_METHOD MATCHES "conan")
    #  Make sure we use GL's own find modules
    list(APPEND CMAKE_MODULE_PATH  ${PROJECT_SOURCE_DIR}/cmake)

    ##################################################################
    ### Install conan-modules/conanfile.txt dependencies          ###
    ### This uses conan to get spdlog,eigen3,h5pp,ceres-solver    ###
    ###    ceres-solver/2.0.0@davidace/development                ###
    ###    h5pp/1.5.1@davidace/stable                             ###
    ###    eigen/3.3.7@davidace/patched                           ###
    ##################################################################




    find_program (
            CONAN_COMMAND
            conan
            HINTS ${CONAN_PREFIX} $ENV{CONAN_PREFIX} ${CONDA_PREFIX} $ENV{CONDA_PREFIX}
            PATHS $ENV{HOME}/anaconda3 $ENV{HOME}/miniconda3 $ENV{HOME}/.conda
            PATH_SUFFIXES bin envs/dmrg/bin
    )
    message(STATUS "Found conan: ${CONAN_COMMAND}")

    # Download cmake-conan automatically, you can also just copy the conan.cmake file
    if(NOT EXISTS "${CMAKE_BINARY_DIR}/conan.cmake")
        message(STATUS "Downloading conan.cmake from https://github.com/conan-io/cmake-conan")
        file(DOWNLOAD "https://github.com/conan-io/cmake-conan/raw/v0.15/conan.cmake"
                "${CMAKE_BINARY_DIR}/conan.cmake")
    endif()

    include(${CMAKE_BINARY_DIR}/conan.cmake)
    conan_add_remote(NAME conan-center       URL https://conan.bintray.com)
    conan_add_remote(NAME conan-community    URL https://api.bintray.com/conan/conan-community/conan)
    conan_add_remote(NAME bincrafters        URL https://api.bintray.com/conan/bincrafters/public-conan)
    conan_add_remote(NAME conan-public INDEX 1 URL https://api.bintray.com/conan/davidace/conan-public)

    if(CMAKE_CXX_COMPILER_ID MATCHES "AppleClang")
        # Let it autodetect libcxx
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
        # There is no libcxx
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        list(APPEND GL_CONAN_SETTINGS SETTINGS compiler.libcxx=libstdc++11)
    endif()
    conan_cmake_run(
            CONANFILE conanfile.txt
            CONAN_COMMAND ${CONAN_COMMAND}
            BUILD_TYPE ${CMAKE_BUILD_TYPE}
            BASIC_SETUP CMAKE_TARGETS
            SETTINGS compiler.cppstd=17
            ${GL_CONAN_SETTINGS}
            ${GL_CONAN_OPTIONS}
            BUILD missing
    )


endif()
