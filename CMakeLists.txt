cmake_minimum_required(VERSION 3.16)
project(OpenCorr LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(OpenCV REQUIRED)
find_package(OpenMP REQUIRED)
find_package(Eigen3 REQUIRED)

file(GLOB OPENCORR_SOURCES "src/*.cpp")
add_library(opencorr STATIC ${OPENCORR_SOURCES})

# OpenCorr's headers use `#include <Eigen>` (the umbrella header with no
# subpath) rather than the conventional `#include <Eigen/Dense>`, so the
# Eigen/ subdirectory itself must be on the include path, not just its
# eigen3/ parent that find_package(Eigen3) normally provides.
get_target_property(EIGEN3_INCLUDE_DIRS Eigen3::Eigen INTERFACE_INCLUDE_DIRECTORIES)

target_include_directories(opencorr PUBLIC
    src
    ${OpenCV_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIRS}/Eigen
    ${CMAKE_CURRENT_SOURCE_DIR}/deps/nanoflann
)

target_link_libraries(opencorr PUBLIC
    ${OpenCV_LIBS}
    Eigen3::Eigen
    OpenMP::OpenMP_CXX
    fftw3
    fftw3f
    fftw3l
)

option(OPENCORR_BUILD_SMOKE_TEST "Build local build-verification smoke tests (*_smoke_test.cpp)" OFF)
if(OPENCORR_BUILD_SMOKE_TEST)
    file(GLOB OPENCORR_SMOKE_TESTS "*_smoke_test.cpp")
    foreach(test_src ${OPENCORR_SMOKE_TESTS})
        get_filename_component(test_name ${test_src} NAME_WE)
        add_executable(${test_name} ${test_src})
        target_link_libraries(${test_name} PRIVATE opencorr)
    endforeach()
endif()
