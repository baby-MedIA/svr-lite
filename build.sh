#!/bin/bash

# Copyright © King's College London
# SPDX-License-Identifier: Apache-2.0

# Set exit on error
set -e

# Set directory locations
sourceDir=$PWD
buildDir=$PWD/build/$mode

# Create build directory if it does not exist
if [ ! -d $buildDir ]; then
  mkdir -p $buildDir
fi

# Build MIRTK if not already present
mirtkSourceDir="$sourceDir/external/mirtk"
mirtkBuildDir="$buildDir/external/mirtk"
if [ ! -d $mirtkBuildDir ]; then
  echo "MIRTK build doesn't exist and must be built"
  mkdir -p $mirtkBuildDir
  echo "Building MIRTK"
  cd $mirtkBuildDir
  cmake -D WITH_TBB="ON" $mirtkSourceDir
  make -j 
fi

# Build SVR-lite
echo "Building SVR-lite"
cd $buildDir
cmake -D WITH_TBB="ON" -D MODULE_SVR="ON" -D CMAKE_BUILD_TYPE=Release -D DEPENDS_MIRTK_DIR="$mirtkBuildDir" ../..

#build
make -j
 