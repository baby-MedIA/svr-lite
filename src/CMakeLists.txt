# Copyright 2008-2017 Imperial College London
# Copyright 2018 Andreas Schuh
# Copyright 2018-2024 King's College London
# SPDX-License-Identifier: Apache-2.0
#
# This file has been modified from the original source.

##############################################################################
# @file  CMakeLists.txt
# @brief Build configuration of MIRTK SVR library.
##############################################################################

set(HEADERS
  ${BINARY_INCLUDE_DIR}/mirtk/SVRExport.h
  ../svr/Application.h
  ../svr/Evaluation.h
  ../svr/ReconstructionParameters.h
  ../svr/ReconstructionModel.h
  ../svr/ReconstructionAlgorithm.h
  ../svr/MeanShift.h
  ../svr/Parallel.h
  ../svr/Utility.h
  ../svr/json.hpp
)

set(SOURCES
  ReconstructionParameters.cc
  ReconstructionModel.cc
  ReconstructionAlgorithm.cc
  Evaluation.cc
  MeanShift.cc
  Utility.cc
)

set(DEPENDS
  LibCommon
  LibNumerics
  LibImage
  LibRegistration
  LibIO
  LibTransformation
  ${Boost_LIBRARIES}
  ${TBB}
)

mirtk_add_library()
