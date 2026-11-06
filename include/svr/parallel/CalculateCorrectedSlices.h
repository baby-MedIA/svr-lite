// Copyright © King's College London
// SPDX-License-Identifier: Apache-2.0
//
// This file has been modified from the original source.

#pragma once

// SVR
#include "svr/ReconstructionModel.h"

using namespace mirtk;
using namespace svr;
using namespace svr::Utility;

namespace svr::Parallel {
    class CalculateCorrectedSlices {
        ReconstructionModel *reconstructor;

    public:
        CalculateCorrectedSlices(ReconstructionModel *reconstructor) : reconstructor(reconstructor) {}

        void operator()(const blocked_range<size_t>& r) const {
            for (size_t inputIndex = r.begin(); inputIndex != r.end(); inputIndex++) {
                //read the current slice
                RealImage slice = reconstructor->_slices[inputIndex];

                //calculate corrected voxels

                for (int t = 0; t < slice.GetT(); t++)
                    for (int j = 0; j < slice.GetY(); j++)
                        for (int i = 0; i < slice.GetX(); i++)
                            if (reconstructor->InSliceMask(inputIndex,i, j))
                                //bias correct and scale the voxel
                                slice(i, j, 0, t) *= exp(-reconstructor->_bias[inputIndex](i, j, 0)) 
                                                        * reconstructor->_scale[inputIndex];
                reconstructor->_corrected_slices[inputIndex] = move(slice);
            } //end of loop for a slice inputIndex
        }

        void operator()() const {
            parallel_for(blocked_range<size_t>(0, reconstructor->_slices.size()), *this);
        }
    };
}