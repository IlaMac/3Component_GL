

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
option(BUILD_SHARED_LIBS     "Enables shared library linking (recommended for MPI)"     ON)
option(GL_ENABLE_TESTS       "Enables CTest testing"                                    ON)
option(GL_ENABLE_H5PP        "h5pp is a C++ front-end to HDF5"                          ON)
option(GL_ENABLE_EIGEN3      "Eigen3 is a linear algebra template library"              ON)
option(GL_ENABLE_SPDLOG      "Spdlog is a fast logging library"                         ON)
option(GL_ENABLE_OPENMP      "Enable OpenMP flags such as -fopenmp"                     ON)
option(GL_ENABLE_MPI         "Enable MPI parallelization"                               ON)
option(GL_ENABLE_LTO         "Enable link time optimization"                            OFF)
option(GL_AS_SUBPROJECT      "Try using this as a subproject with add_subdirectory()"   OFF)
option(GL_DEPS_IN_SUBDIR     "Install dependencies into CMAKE_INSTALL_PREFIX/<libname>" ON)
option(GL_PRINT_INFO         "Prints additional information during CMake configure"     ON)

if(NOT GL_ENABLE_MPI)
    message(FATAL_ERROR "This project uses MPI and should be enabled by default.")
endif()

if(GL_ENABLE_MPI AND NOT BUILD_SHARED_LIBS)
    message(WARNING "This project uses MPI which strongly encourages shared library linking")
endif()
#################################################################
### Make an "enum" for valid download methods                 ###
#################################################################
set(GL_VALID_DOWNLOAD_METHODS none find fetch find-or-fetch conan)
set(GL_DOWNLOAD_METHOD fetch CACHE STRING "Download method for external dependencies")
set_property(CACHE GL_DOWNLOAD_METHOD  PROPERTY STRINGS ${VALID_DOWNLOAD_METHODS})
if(NOT GL_DOWNLOAD_METHOD IN_LIST GL_VALID_DOWNLOAD_METHODS)
    message(FATAL_ERROR "GL_DOWNLOAD_METHOD must be one of ${GL_VALID_DOWNLOAD_METHODS}")
endif()


# Used when GL is included as subproject (e.g., as Git submodule/subtree) in the source
# tree of a project that uses it. Users may set the non-cached variable GL_AS_SUBPROJECT
# before add_subdirectory(GL)
if(NOT GL_AS_SUBPROJECT)
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


############################################################
### Print extra info during CMake configure              ###
############################################################
include(cmake/PrintHostInfo.cmake)
include(cmake/PrintOptionInfo.cmake)



### Print compiler details
message("Compiling with: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")


### Create a helper target with all the options you need to compile
### This one is useful to pass to tests and examples so that they compile
### In the same way as our main project.
add_library(project-settings INTERFACE)
target_compile_features(project-settings INTERFACE cxx_std_17)



# Set default install directory for dependencies
# Append directory <libname> to CMAKE_INSTALL_PREFIX
# Useful if you want to have separate directories for each libs
# and to easily delete them individually
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT AND NOT GL_IS_SUBPROJECT)
    message(STATUS "Setting default install prefix -- CMAKE_INSTALL_PREFIX --> ${CMAKE_BINARY_DIR}/deps-install")
    set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/deps-install" CACHE PATH "..." FORCE)
endif ()



### Add all source files recursively
file(GLOB_RECURSE SOURCES "source/*.cpp")
file(GLOB_RECURSE HEADERS "source/*.h")
add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS} source/class_tic_toc.cpp source/class_tic_toc.h)



##################################################################
###  Force cmake to find .a library suffixes on static builds  ###
##################################################################
if(BUILD_SHARED_LIBS)
    set(LINK_TYPE SHARED)
else()
    set(LINK_TYPE STATIC)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_STATIC_LIBRARY_SUFFIX} )
    #target_link_libraries  (project-settings INTERFACE -static)   ### Static linkage
endif()


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





# Fetch or find all the dependencies
include(cmake/SetupDependenciesFind.cmake)
include(cmake/SetupDependenciesFetch.cmake)
include(cmake/SetupDependenciesConan.cmake)

##################################################################
### Link all the things!                                       ###
##################################################################
target_link_libraries(project-settings INTERFACE  ${CONAN_TARGETS} ${NATIVE_TARGETS} Threads::Threads) # Order is important here
target_link_libraries(${PROJECT_NAME} PRIVATE project-settings)  # Link everything to main project




# Print summary of CMake configuration
if (GL_PRINT_INFO)
    include(cmake/PrintTargetInfo.cmake)
    include(cmake/getExpandedTarget.cmake)
    expand_target_all_targets(project-settings GL_TARGETS_EXPANDED)
    message(STATUS "| GL TARGET SUMMARY")
    message(STATUS "|--------------------")
    foreach(tgt ${GL_TARGETS_EXPANDED})
        print_target_info(${tgt} "| ")
    endforeach()
    print_target_info(${PROJECT_NAME} "| ")
endif()




# Enable link time optimization
if(GL_ENABLE_LTO)
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



############################################################
### Set  the same microarchitecture for c++ and OpenBLAS ###
############################################################


if(GL_MARCH)
    if (GL_MARCH MATCHES "zen")
        set(CXX_MARCH znver1)
    else()
        string(TOLOWER ${GL_MARCH} CXX_MARCH)
    endif()
    message(STATUS "CPU optimizations: ${CXX_MARCH}")
else()
    set(GL_MARCH "native")
    set(CXX_MARCH "native")
endif()


###########################################
###  Apply RELEASE/DEBUG compile flags  ###
###########################################
# Segfault from -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC when using ceres, but -D_GLIBCXX_DEBUG fixes gdb incomplete types
target_compile_options(project-settings INTERFACE $<$<COMPILE_LANGUAGE:CXX>: -Wall -march=${CXX_MARCH} -mtune=${CXX_MARCH}>)
target_compile_options(project-settings INTERFACE $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:RELEASE>>:          -O3 -g -DNDEBUG>)
target_compile_options(project-settings INTERFACE $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:DEBUG>>:            -O0 -g3 -Wpedantic -Wextra -fstack-protector -D_FORTIFY_SOURCE=2 -D_GLIBCXX_DEBUG>)
target_compile_options(project-settings INTERFACE $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:RELWITHDEBINFO>>:   -O1 -g3 -Wpedantic -Wextra -fstack-protector -D_FORTIFY_SOURCE=2>)
target_compile_options(project-settings INTERFACE $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:PROFILE>>:          -O3 -g  -DNDEBUG  -lprofiler  -ftime-report>)



### CTest testing
if(GL_ENABLE_TESTS AND NOT GL_AS_SUBPROJECT)
    enable_testing()
    add_subdirectory(tests)
endif()


