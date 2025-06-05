cmake_minimum_required(VERSION 3.10)
project(geffe_generator)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(SOURCE_DIR ${CMAKE_SOURCE_DIR}/src)
set(MAIN_DIR ${CMAKE_SOURCE_DIR}/main)
set(TEST_DIR ${CMAKE_SOURCE_DIR}/test)

set(SOURCES
    ${SOURCE_DIR}/geffe_generator.cpp
    ${SOURCE_DIR}/lfsr.cpp
    ${SOURCE_DIR}/register_recovery.cpp
)

include_directories(${SOURCE_DIR})

add_executable(geffe_generator ${MAIN_DIR}/geffe_generator_driver.cpp ${SOURCES})
add_executable(test_L1 ${TEST_DIR}/L1_recovery_test.cpp ${SOURCES})
add_executable(test_L2 ${TEST_DIR}/L2_recovery_test.cpp ${SOURCES})