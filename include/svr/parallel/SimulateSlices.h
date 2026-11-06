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

    class SimulateSlices {
        ReconstructionModel *reconstructor;

    public:
        SimulateSlices(ReconstructionModel *reconstructor) : reconstructor(reconstructor) {}

        void operator()(const blocked_range<size_t>& r) const {
            for (size_t inputIndex = r.begin(); inputIndex != r.end(); inputIndex++) {
                const RealImage& slice = reconstructor->_slices[inputIndex];
                //Calculate simulated slice
                reconstructor->_simulated_slices[inputIndex].Initialize(slice.Attributes());
                reconstructor->_simulated_weights[inputIndex].Initialize(slice.Attributes());
                reconstructor->_simulated_inside[inputIndex].Initialize(slice.Attributes());
                reconstructor->_slice_inside[inputIndex] = false;

                for (int l = 0; l < reconstructor->_reconstructed.GetT(); l++) 
                    for (size_t i = 0; i < reconstructor->_volcoeffs[inputIndex].size(); i++) // X
                        for (size_t j = 0; j < reconstructor->_volcoeffs[inputIndex][i].size(); j++) {// Y
                            if (reconstructor->InSliceMask(inputIndex, i, j) || reconstructor->_no_masking_background) { 
                                double weight = 0;
                                for (size_t k = 0; k < reconstructor->_volcoeffs[inputIndex][i][j].size(); k++) { // coeffs size
                                    const POINT3D& p = reconstructor->_volcoeffs[inputIndex][i][j][k];
                                    reconstructor->_simulated_slices[inputIndex](i, j, 0, l) += p.value * reconstructor->_reconstructed(p.x, p.y, p.z, l);
                                    weight += p.value;
                                    if (reconstructor->_reconstructed_mask(p.x, p.y, p.z) == 1) {
                                        reconstructor->_simulated_inside[inputIndex](i, j, 0, l) = 1;
                                        reconstructor->_slice_inside[inputIndex] = true;
                                    }
                                }
                                if (weight > 0) {
                                    reconstructor->_simulated_slices[inputIndex](i, j, 0, l) /= weight;
                                    reconstructor->_simulated_weights[inputIndex](i, j, 0, l) = weight;
                                }
                            }
                        }
            }
        }

        void operator()() const {
            parallel_for(blocked_range<size_t>(0, reconstructor->_slices.size()), *this);
        }
    };

}