add_executable(rtcqt-resampling-tests
    resampling_tests.cpp
)

target_include_directories(rtcqt-resampling-tests
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
)

target_compile_features(rtcqt-resampling-tests PRIVATE cxx_std_17)

add_test(
    NAME rtcqt-resampling-tests
    COMMAND rtcqt-resampling-tests
)

set(RTCQT_PFFFT_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/../submodules/pffft/pffft.c
    ${CMAKE_CURRENT_SOURCE_DIR}/../submodules/pffft/pffft_common.c
    ${CMAKE_CURRENT_SOURCE_DIR}/../submodules/pffft/pffft_double.c
)

add_executable(rtcqt-cqt-integration-tests
    cqt_integration_tests.cpp
    ${RTCQT_PFFFT_SOURCES}
)

target_include_directories(rtcqt-cqt-integration-tests
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
)

target_compile_features(rtcqt-cqt-integration-tests PRIVATE cxx_std_17)

add_test(
    NAME rtcqt-cqt-integration-tests
    COMMAND rtcqt-cqt-integration-tests
)
