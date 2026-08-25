cmake_minimum_required(VERSION 3.15)
project(taplite_calibration_native LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(TAPLITE_CALIBRATION_SOURCES
    src/AutoCalibration.cpp)

find_package(OpenMP)

option(TAPLITE_BUILD_PYTHON "Build the pybind11 extension" ON)
if(TAPLITE_BUILD_PYTHON)
    find_package(pybind11 CONFIG REQUIRED)
    pybind11_add_module(_native
        binding.cpp
        src/TAPLite.cpp
        ${TAPLITE_CALIBRATION_SOURCES})
    target_include_directories(_native PRIVATE src)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(_native PRIVATE OpenMP::OpenMP_CXX)
    endif()
endif()

option(TAPLITE_BUILD_CPP_TESTS "Build focused C++ calibration tests" OFF)
if(TAPLITE_BUILD_CPP_TESTS)
    enable_testing()
    foreach(test_name IN ITEMS test_auto_calibration)
        add_executable(${test_name} tests/${test_name}.cpp ${TAPLITE_CALIBRATION_SOURCES})
        target_include_directories(${test_name} PRIVATE src)
        if(OpenMP_CXX_FOUND)
            target_link_libraries(${test_name} PRIVATE OpenMP::OpenMP_CXX)
        endif()
        add_test(NAME ${test_name} COMMAND ${test_name})
    endforeach()
endif()
