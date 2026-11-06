# Copyright 2008-2017 Imperial College London
# Copyright 2018-2024 King's College London
# SPDX-License-Identifier: Apache-2.0
#
# This file has been modified from the original source.

cmake_minimum_required(VERSION 2.8.12 FATAL_ERROR)
project("SVR-lite")
set(CMAKE_CXX_STANDARD 17)

## Dependencies ##

list(INSERT CMAKE_MODULE_PATH 0 "${CMAKE_CURRENT_SOURCE_DIR}/config")
if (NOT COMMAND mirtk_configure_module)
  find_package(MIRTK REQUIRED COMPONENTS CMake)
endif()

find_package(spdlog REQUIRED)
find_package(OpenMP)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

mirtk_configure_module()

#add_subdirectory(test)