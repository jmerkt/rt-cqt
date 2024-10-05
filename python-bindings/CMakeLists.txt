cmake_minimum_required(VERSION 3.4...3.18)
project(rtcqt-python)

add_subdirectory(../submodules/pybind11 ${CMAKE_CURRENT_BINARY_DIR}/pybind11)
pybind11_add_module(rtcqt-python src/main.cpp)

# EXAMPLE_VERSION_INFO is defined by setup.py and passed into the C++ code as a
# define (VERSION_INFO) here.
set(EXAMPLE_VERSION_INFO 0.0.1)
target_compile_definitions(rtcqt-python PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})

# rt-cqt source files
set(ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/..)

set(SOURCE_FILES
	${ROOT_DIR}/submodules/pffft/pffft.c
	${ROOT_DIR}/submodules/pffft/pffft_common.c
	${ROOT_DIR}/submodules/pffft/pffft_double.c
)

target_sources(rtcqt-python PRIVATE ${SOURCE_FILES})

add_compile_definitions(INCLUDE_PYTHON_BINDING=true)
