cmake_minimum_required(VERSION 3.10)
project(PrivLR)
add_subdirectory(src)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

option(PrivLR_TEST "Build tests" ON)
if(PrivLR_TEST)
  message(STATUS "Tests enabled")
  add_subdirectory(test)
endif()

set(INSTALL_DIR "${CMAKE_INSTALL_PREFIX}")
set(BUILD_DIR "${PROJECT_SOURCE_DIR}/build")
