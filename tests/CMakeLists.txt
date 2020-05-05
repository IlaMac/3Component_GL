cmake_minimum_required(VERSION 3.12)
project(tests)
add_custom_target(all-tests)

if(TARGET Eigen3::Eigen OR TARGET CONAN_PKG::Eigen3)
    list(APPEND TESTS ${PROJECT_SOURCE_DIR}/eigen.cpp)
endif()
if(TARGET h5pp::h5pp OR TARGET CONAN_PKG::h5pp)
    list(APPEND TESTS ${PROJECT_SOURCE_DIR}/h5pp.cpp)
endif()
if(TARGET openmp::openmp)
    list(APPEND TESTS ${PROJECT_SOURCE_DIR}/openmp.cpp)
endif()

if(TESTS)
    foreach (test ${TESTS})
        get_filename_component(test_src ${test} NAME)
        get_filename_component(test_nwe ${test} NAME_WE)
        add_executable(${test_nwe} ${test_src} ${test_inc})
        target_link_libraries(${test_nwe} PRIVATE project-settings)
        add_test(NAME ${test_nwe} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND ${test_nwe})
        if (MSVC)
            set_target_properties(${test_nwe} PROPERTIES LINK_FLAGS "/ignore:4099")
        endif()
        add_dependencies(all-tests ${test_nwe})
    endforeach ()
    add_custom_command(
            TARGET all-tests
            POST_BUILD
            COMMENT "Running Tests"
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            DEPENDS project-settings
            COMMAND ${CMAKE_CTEST_COMMAND} -C $<CONFIG> --output-on-failure)
endif()