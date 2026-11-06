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
        /// Class for slice scale calculation
    class Scale {
        ReconstructionModel *reconstructor;

    public:
        Scale(ReconstructionModel *reconstructor) : reconstructor(reconstructor) {}

        void operator()(const blocked_range<size_t>& r) const {
            for (size_t inputIndex = r.begin(); inputIndex != r.end(); inputIndex++) {
                //alias the current bias image
                const RealImage& b = reconstructor->_bias[inputIndex];

                //initialise calculation of scale
                double scalenum = 0;
                double scaleden = 0;

                //voxelfunction candidate
                for (int l = 0; l < reconstructor->_slices[inputIndex].Attributes().T(); l++ )
                    for (int j = 0; j < reconstructor->_slices[inputIndex].GetY(); j++)
                        for (int i = 0; i < reconstructor->_slices[inputIndex].GetX(); i++)
                            if (reconstructor->InSliceMask(inputIndex,i, j)) { 
                                if (reconstructor->_simulated_weights[inputIndex](i, j, 0, l) > 0.99) {
                                    //scale - intensity matching
                                    const double eb = exp(-b(i, j, 0));
                                    scalenum += reconstructor->_weights[inputIndex](i, j, 0, l) 
                                                * reconstructor->_slices[inputIndex](i, j, 0, l) 
                                                * eb * reconstructor->_simulated_slices[inputIndex](i, j, 0, l);
                                    scaleden += reconstructor->_weights[inputIndex](i, j, 0, l) 
                                                * reconstructor->_slices[inputIndex](i, j, 0, l) 
                                                * eb * reconstructor->_slices[inputIndex](i, j, 0, l) * eb;
                                }
                            }

                //calculate scale for this slice
                reconstructor->_scale[inputIndex] = scaleden > 0 ? scalenum / scaleden : 1;
            } //end of loop for a slice inputIndex
        }

        void operator()() const {
            parallel_for(blocked_range<size_t>(0, reconstructor->_slices.size()), *this);
        }
    };
}