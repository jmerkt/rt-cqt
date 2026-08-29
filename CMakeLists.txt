cmake_minimum_required(VERSION 3.4...3.18)
project(rtcqt)

option(RTCQT_BUILD_EXAMPLES "Build the rt-cqt examples" ON)
option(RTCQT_BUILD_PYTHON_BINDINGS "Build the Python bindings" ON)

if(RTCQT_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

if(RTCQT_BUILD_PYTHON_BINDINGS)
    add_subdirectory(python-bindings)
endif()

include(CTest)
if(BUILD_TESTING)
    add_subdirectory(tests)
endif()
