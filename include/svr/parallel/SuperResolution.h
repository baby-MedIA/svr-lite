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
    class Superresolution {
        ReconstructionModel *reconstructor;

    public:
        RealImage confidence_map;
        RealImage addon;

        Superresolution(ReconstructionModel *reconstructor) : reconstructor(reconstructor) {
            //Clear addon
            addon.Initialize(reconstructor->_reconstructed.Attributes());

            //Clear confidence map
            confidence_map.Initialize(reconstructor->_reconstructed.Attributes());
        }

        Superresolution(Superresolution& x, split) : Superresolution(x.reconstructor) {}

        void operator()(const blocked_range<size_t>& r) {
            RealImage slice;

            // For each slice
            for (size_t sliceIndex = r.begin(); sliceIndex < r.end(); sliceIndex++) {
                if (reconstructor->_volcoeffs[sliceIndex].empty())
                    continue;

                // read the current slice and mask
                slice = reconstructor->_slices[sliceIndex];

                //Update reconstructed volume using current slice
                //Distribute error to the volume
                for (int l = 0; l < reconstructor->_reconstructed.GetT(); l++) {
                    for (size_t i = 0; i < reconstructor->_volcoeffs[sliceIndex].size(); i++) {
                        for (size_t j = 0; j < reconstructor->_volcoeffs[sliceIndex][i].size(); j++) {
                            // If in ROI, apply super-resolution functionality
                            if (reconstructor->InSliceMask(sliceIndex,i, j)) {
                            
                                //bias correct and scale the slice
                                slice(i, j, 0, l) *= exp(-reconstructor->_bias[sliceIndex](i, j, 0)) * reconstructor->_scale[sliceIndex];

                                if (reconstructor->_simulated_slices[sliceIndex](i, j, 0, l) > 0)
                                    slice(i, j, 0, l) -= reconstructor->_simulated_slices[sliceIndex](i, j, 0, l);
                                else
                                    slice(i, j, 0, l) = 0;

                                #pragma omp simd
                                for (size_t k = 0; k < reconstructor->_volcoeffs[sliceIndex][i][j].size(); k++) {
                                    const POINT3D& p = reconstructor->_volcoeffs[sliceIndex][i][j][k];
                                    const auto multiplier = reconstructor->_robust_slices_only ? 1 : reconstructor->_weights[sliceIndex](i, j, 0, l);

                                    double ssim_weight = 1;
                                    if (reconstructor->_structural)
                                        ssim_weight = reconstructor->_slice_ssim_maps[sliceIndex](i, j, 0, l);

                                    addon(p.x, p.y, p.z, l) += p.value * slice(i, j, 0, l) * ssim_weight * multiplier * reconstructor->_slice_weight[sliceIndex];
                                    confidence_map(p.x, p.y, p.z, l) += p.value * ssim_weight * multiplier * reconstructor->_slice_weight[sliceIndex];
                                }
                            }
                        }
                    }
                }
            }
        }

        void join(const Superresolution& y) {
            addon += y.addon;
            confidence_map += y.confidence_map;
        }

        void operator()() {
            parallel_reduce(blocked_range<size_t>(0, reconstructor->_slices.size()), *this);
        }
    };

}