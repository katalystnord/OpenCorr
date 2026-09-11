cmake_minimum_required(VERSION 3.16)
project(OpenCorr LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(OpenCV REQUIRED)
find_package(OpenMP REQUIRED)
find_package(Eigen3 REQUIRED)

file(GLOB OPENCORR_SOURCES "src/*.cpp")
# Coverage instrumentation, for tools/coverage.sh.
#
# ⚑ This block must stay ABOVE add_library() below. add_compile_options() only
# affects targets created after it, and when this sat further down the library
# was built uninstrumented while only the smoke tests were -- so gcovr reported
# "100% of 83 lines", which is the headers and nothing else. A coverage number
# that is impossibly high is a measurement fault, not a result. Off by default and built into
# a separate directory: --coverage changes the object code, so measuring what the
# tests reach and producing the library anyone links are different jobs.
option(OPENCORR_COVERAGE "Build with gcov instrumentation (see tools/coverage.sh)" OFF)
if(OPENCORR_COVERAGE)
    # -O0 so line and branch counts correspond to the source as written; an
    # optimised build merges and elides lines, and the report then describes
    # something the reader cannot point at.
    add_compile_options(--coverage -O0 -g)
    add_link_options(--coverage)
endif()

# ⚑ ADDRESS AND UNDEFINED-BEHAVIOUR SANITIZERS. This is the tool the mutation
# work in SurView and PatternFab pointed at: of the mutants left standing there
# after every case worth writing had been written, the large majority read or
# write ONE ELEMENT past the end of a vector, and no assertion can catch that --
# the value read is whatever sits beyond the buffer and nothing downstream need
# change. A sanitizer catches the access itself rather than its consequences,
# which is the only way that class gets closed.
#
# It matters more here than there. This library is index and pointer arithmetic
# in OpenMP hot loops, and a defect in it is a WRONG NUMBER rather than a wrong
# message. Three crash bugs in hypercine's read_header() were found by reading
# the code (fork issue #18); a sanitizer finds that shape without anyone having
# to suspect it first.
#
# -fno-sanitize-recover is what makes it a test rather than a log: without it
# UBSan prints a diagnostic, carries on, and the suite stays green while telling
# the truth in a place nobody reads. Off by default and built into its own
# directory, for the same reason coverage is: the instrumentation changes the
# object code, so producing the library anyone links and checking it are
# different jobs.
option(OPENCORR_SANITIZE "Build with address and undefined-behaviour sanitizers" OFF)
if(OPENCORR_SANITIZE)
    add_compile_options(-fsanitize=address,undefined -fno-sanitize-recover=all
                        -fno-omit-frame-pointer -g -O1)
    add_link_options(-fsanitize=address,undefined)
endif()

add_library(opencorr STATIC ${OPENCORR_SOURCES} deps/hypercine/hypercine.cpp)

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

# PUBLIC, not PRIVATE: oc_cine.h only forward-declares hypercine::HyperCine so
# normal consumers never need this, but a downstream consumer (SurView) reaching
# for cine metadata Cine2D's minimal API doesn't expose (frame rate, trigger time,
# etc.) may need to include hypercine.h directly -- if USE_FLOAT_STORAGE weren't
# propagated to them too, storage_t would be undefined in their translation unit,
# or worse, silently resolve to a different type than what libopencorr.a was
# actually compiled with (an ODR/ABI mismatch).
target_include_directories(opencorr PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/deps/hypercine
)

# hypercine's storage_t is a compile-time choice (see deps/hypercine/hypercine.h);
# float matches Image2D's Eigen::MatrixXf representation.
target_compile_definitions(opencorr PUBLIC USE_FLOAT_STORAGE=1)

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
    # Registered with CTest, not merely built. Built-but-unregistered tests have
    # to be found and run one at a time by whoever remembers they exist, which
    # is how a suite quietly stops being run in full: `ctest` reported nothing
    # to do here while twenty smoke tests sat beside it.
    enable_testing()

    file(GLOB OPENCORR_SMOKE_TESTS "*_smoke_test.cpp")
    foreach(test_src ${OPENCORR_SMOKE_TESTS})
        get_filename_component(test_name ${test_src} NAME_WE)
        add_executable(${test_name} ${test_src})
        target_link_libraries(${test_name} PRIVATE opencorr)
        add_test(NAME ${test_name} COMMAND ${test_name})
        # Several read example images by a path relative to the repository
        # root, so that is where they are run from regardless of build layout.
        set_tests_properties(${test_name} PROPERTIES
            WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
    endforeach()
endif()
