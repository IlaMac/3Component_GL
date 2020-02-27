

### Enforce an out-of-source build, because c'mon...
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(STATUS "CMakeTemplate requires an out-of-source build.")
    message(STATUS "Please remove these files from ${CMAKE_BINARY_DIR} first:")
    message(STATUS "CMakeCache.txt")
    message(STATUS "CMakeFiles")
    message(STATUS "Once these files are removed, create a separate directory")
    message(STATUS "and run CMake from there")
    message(FATAL_ERROR "in-source build detected")
endif()

cmake_minimum_required(VERSION 3.12)
### Setup project name and language
project(GL_3component CXX)

### Set options for compilation
option(ENABLE_TESTS                 "Enables CTest testing"                                                           OFF)
option(BUILD_SHARED_LIBS            "Sets shared linking on/off"                                                      ON)
option(ENABLE_EIGEN3                "Use Eigen3 Library. Guide with -DEIGEN_INCLUDE_DIR:PATH=<dir-to-eigen3>"         OFF)
option(ENABLE_SPDLOG                "Use SpdLog logging library. Guide with -Dspdlog_DIR:PATH=<dir-to-spdlog>"        OFF)
option(ENABLE_H5PP                  "Use h5pp, a C++ front-end to HDF5"                                               OFF)
option(ENABLE_OPENMP                "Enable OpenMP flags such as -fopenmp"                                            ON)
option(ENABLE_LTO                   "Enable link time optimization"                                                   OFF)
option(AS_SUBPROJECT                "Try using this as a subproject with add_subdirectory()"                          OFF)
option(APPEND_LIBSUFFIX             "Append <libname> to CMAKE_INSTALL_PREFIX"                                        ON)



#################################################################
### Make an "enum" for valid download methods                 ###
#################################################################
set(VALID_DOWNLOAD_METHODS none find-only conan native)
set(DOWNLOAD_METHOD find-only CACHE STRING "Download method for external dependencies")
set_property(CACHE DOWNLOAD_METHOD  PROPERTY STRINGS ${VALID_DOWNLOAD_METHODS})
if(NOT DOWNLOAD_METHOD IN_LIST VALID_DOWNLOAD_METHODS)
    message(FATAL_ERROR "DOWNLOAD_METHOD must be one of ${VALID_DOWNLOAD_METHODS}")
endif()


### h5pp includes eigen and spdlog as a dependencies
if(ENABLE_H5PP)
    message(STATUS "NOTE: Eigen3 and spdlog are included with h5pp")
    set(ENABLE_EIGEN3 OFF)
    set(ENABLE_SPDLOG OFF)
endif()




# Used when h5pp is included as subproject (e.g., as Git submodule/subtree) in the source
# tree of a project that uses it. Users may set the non-cached variable AS_SUBPROJECT
# before add_subdirectory(h5pp)
if(NOT AS_SUBPROJECT)
    if("^${CMAKE_SOURCE_DIR}$" STREQUAL "^${PROJECT_SOURCE_DIR}$")
        set (AS_SUBPROJECT FALSE)
    else()
        message(STATUS "Detected usage as subproject")
        set (AS_SUBPROJECT TRUE)
    endif()
endif()


################################################################
### Get git version number                                   ###
### Generates a header gitversion/gitversion.h               ###
### Include it using #include <gitversion.h>                 ###
### Gives a namespace GIT:: with several git identifiers     ###
################################################################
include(cmake/gitversion.cmake)


### Print operating system details
include(cmake/PrintHostInfo.cmake)

### Print compiler details
message("Compiling with: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")


### Create a helper target with all the options you need to compile
### This one is useful to pass to tests and examples so that they compile
### In the same way as our main project.
add_library(project-settings INTERFACE)
target_compile_features(project-settings INTERFACE cxx_std_17)


# Let cmake find our Find<package>.cmake modules
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# Require the filesystem header (and possibly stdc++fs library)
find_package(Filesystem COMPONENTS Final Experimental)
if (TARGET std::filesystem)
    target_link_libraries(project-settings INTERFACE std::filesystem)
elseif (DOWNLOAD_METHOD MATCHES "find|none")
    message(STATUS "Your compiler lacks std::filesystem. Set DOWNLOAD_METHOD to 'conan' or 'native' to get the ghc::filesystem replacement")
    message(STATUS "Read more about ghc::filesystem here: https://github.com/gulrak/filesystem")
    message(FATAL_ERROR "<filesystem> header and/or library not found")
else ()
    message(STATUS "Your compiler lacks std::filesystem. Getting ghc::filesystem replacement")
    message(STATUS "Read more about ghc::filesystem here: https://github.com/gulrak/filesystem")
    include(cmake/Fetch_ghcFilesystem.cmake)
    target_link_libraries(project-settings INTERFACE ghcFilesystem::ghc_filesystem)
    list(APPEND H5PP_DIRECTORY_HINTS ${ghcfilesystem_install_prefix})
endif ()


### Add all source files recursively
file(GLOB_RECURSE SOURCES "source/*.cpp")
file(GLOB_RECURSE HEADERS "source/*.h")
add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})



##################################################################
### Preempt Threads::Threads                                   ###
### It's looked for in dependencies, so we make it right       ###
### before it's done wrong, i.e. with pthread instead of       ###
### -lpthread.                                     ###
##################################################################
if(NOT TARGET Threads::Threads)
    set(CMAKE_THREAD_PREFER_PTHREAD FALSE)
    set(THREADS_PREFER_PTHREAD_FLAG FALSE)
    find_package(Threads)
    set_target_properties(Threads::Threads PROPERTIES INTERFACE_LINK_LIBRARIES pthread)
endif()


### Link MPI
find_package(MPI REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE MPI::MPI_CXX)




# Set default install directory for dependencies
# Append directory <libname> to CMAKE_INSTALL_PREFIX
# Useful if you want to have separate directories for each libs
# and to easily delete them individually
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT AND NOT H5PP_IS_SUBPROJECT)
    message(STATUS "Setting default install prefix -- CMAKE_INSTALL_PREFIX --> ${CMAKE_BINARY_DIR}/deps-install")
    set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/deps-install" CACHE PATH "..." FORCE)
endif ()
if (APPEND_LIBSUFFIX)
    set(H5PP_INSTALL_DIR h5pp/)
    set(H5PP_CONFIG_DIR h5pp/share/h5pp/cmake)
    set(eigen3_install_prefix ${CMAKE_INSTALL_PREFIX}/Eigen3)
    set(spdlog_install_prefix ${CMAKE_INSTALL_PREFIX}/spdlog)
    set(hdf5_install_prefix ${CMAKE_INSTALL_PREFIX}/hdf5)
    set(ghcfilesystem_install_prefix ${CMAKE_INSTALL_PREFIX}/ghcFilesystem)
else ()
    set(H5PP_CONFIG_DIR share/h5pp/cmake)
    set(eigen3_install_prefix ${CMAKE_INSTALL_PREFIX})
    set(spdlog_install_prefix ${CMAKE_INSTALL_PREFIX})
    set(hdf5_install_prefix ${CMAKE_INSTALL_PREFIX})
    set(ghcfilesystem_install_prefix ${CMAKE_INSTALL_PREFIX})
endif ()




# Fetch or find all the dependencies
if(DOWNLOAD_METHOD MATCHES "conan")
    include(cmake/SetupDependenciesConan.cmake)
elseif(DOWNLOAD_METHOD MATCHES "find|native")
    include(cmake/SetupDependenciesNative.cmake)
endif()
include(cmake/PrintTargetInfo.cmake)
target_link_libraries(${PROJECT_NAME} PRIVATE project-settings)  # Link everything to main project
print_target_info(${PROJECT_NAME})


# Enable link time optimization
if(ENABLE_LTO)
    ### Check if we can enable LTO
    include(CheckIPOSupported)
    check_ipo_supported(RESULT supported OUTPUT error)
    if(supported)
        message(STATUS "IPO / LTO supported, enabling")
        set_target_properties(${PROJECT_NAME} PROPERTIES INTERPROCEDURAL_OPTIMIZATION ON)
    else()
        message(STATUS "IPO / LTO not supported: <${error}>, will not enable")
    endif()
endif()



### CTest testing
if(ENABLE_TESTS AND NOT AS_SUBPROJECT)
    enable_testing()
    add_custom_target(all-tests)
    if(TARGET h5pp::h5pp OR TARGET CONAN_PKG::h5pp)
        add_subdirectory(tests/h5pp)
        add_dependencies(all-tests test-h5pp)
    endif()
    if(TARGET spdlog::spdlog OR TARGET CONAN_PKG::spdlog)
        add_subdirectory(tests/spdlog)
        add_dependencies(all-tests test-spdlog)
    endif()
    if(TARGET Eigen3::Eigen OR TARGET CONAN_PKG::Eigen3 OR TARGET CONAN_PKG::eigen)
        add_subdirectory(tests/eigen3)
        add_dependencies(all-tests test-eigen3)
    endif()
    if(TARGET openmp::openmp)
        add_subdirectory(tests/openmp)
        add_dependencies(all-tests test-openmp)
    endif()
    #Run all tests as soon as the tests have been built
    add_custom_command(TARGET all-tests
            POST_BUILD
            COMMENT "Running Tests"
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
            COMMAND ${CMAKE_CTEST_COMMAND} -C $<CONFIGURATION> --output-on-failures)

endif()


