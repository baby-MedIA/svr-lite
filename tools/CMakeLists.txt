# Copyright 2008-2017 Imperial College London
# Copyright 2018-2024 King's College London
# SPDX-License-Identifier: Apache-2.0
#
# This file has been modified from the original source.

##############################################################################
# @file  CMakeLists.txt
# @brief Build configuration of MIRTK SVR tools.
##############################################################################

mirtk_add_executable(
  svr
  SOURCES
    svr.cc
  DEPENDS
    LibCommon
    LibNumerics
    LibImage
    LibIO
    LibRegistration
    LibTransformation
    LibSVR
    ${TBB}
)

mirtk_add_executable(
  reconstruct
  SOURCES
    reconstruct.cc
  DEPENDS
    LibCommon
    LibNumerics
    LibImage
    LibIO
    LibRegistration
    LibTransformation
    LibSVR
    ${TBB}
)

mirtk_add_executable(
  evaluate-slice
  SOURCES
    evaluate-slice.cc
  DEPENDS
    LibCommon
    LibNumerics
    LibImage
    LibIO
    LibRegistration
    LibTransformation
    LibSVR
    ${TBB}
)


mirtk_add_executable(
  ndiff
  SOURCES
    ndiff.cc
  DEPENDS
    LibCommon
    LibNumerics
    LibImage
    LibIO
    LibRegistration
    LibTransformation
    LibSVR
    ${TBB}
)