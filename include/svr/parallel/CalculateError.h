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

    class CalculateError {
        ReconstructionModel *reconstructor;

    public:
        CalculateError(ReconstructionModel *reconstructor) : reconstructor(reconstructor) {}

        void operator()(const blocked_range<size_t>& r) const {
            RealImage slice;

            for (size_t inputIndex = r.begin(); inputIndex != r.end(); inputIndex++) {
                //read the current slice
                slice = reconstructor->_slices[inputIndex];

                //initialise
                reconstructor->_error[inputIndex].Initialize(slice.Attributes());

                //calculate error
                for (int t = 0; t < slice.GetT(); t++)
                    for (int j = 0; j < slice.GetY(); j++)
                        for (int i = 0; i < slice.GetX(); i++)
                        if (reconstructor->InSliceMask(inputIndex,i, j))
                            if (reconstructor->_simulated_weights[inputIndex](i, j, 0, t) > 0) {
                                //bias correct and scale the voxel
                                slice(i, j, 0, t) *= exp(-reconstructor->_bias[inputIndex](i, j, 0)) 
                                                    * reconstructor->_scale[inputIndex];
                                //subtract simulated voxel
                                slice(i, j, 0, t) -= reconstructor->_simulated_slices[inputIndex](i, j, 0, t);
                                //assign as error
                                reconstructor->_error[inputIndex](i, j, 0, t) = slice(i, j, 0, t);
                            }

            } //end of loop for a slice inputIndex
        }

        void operator()() const {
            parallel_for(blocked_range<size_t>(0, reconstructor->_slices.size()), *this);
        }
    };
}