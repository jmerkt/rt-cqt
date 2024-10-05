cmake_minimum_required(VERSION 3.4...3.18)
project(rtcqt-examples)

# rt-cqt source files
set(ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/..)

set(PFFFT_SOURCE_FILES
	${ROOT_DIR}/submodules/pffft/pffft.c
	${ROOT_DIR}/submodules/pffft/pffft_common.c
	${ROOT_DIR}/submodules/pffft/pffft_double.c
)

add_executable(${PROJECT_NAME}-cqt
	cqt.cpp
	${PFFFT_SOURCE_FILES}
)

target_include_directories(${PROJECT_NAME}-cqt
	PUBLIC
	$<BUILD_INTERFACE:${ROOT_DIR}/include>
)

target_link_libraries(${PROJECT_NAME}-cqt
	PRIVATE
)

add_executable(${PROJECT_NAME}-scqt
	scqt.cpp
)

target_include_directories(${PROJECT_NAME}-scqt
	PUBLIC
	$<BUILD_INTERFACE:${ROOT_DIR}/include>
)

target_link_libraries(${PROJECT_NAME}-scqt
	PRIVATE
)
