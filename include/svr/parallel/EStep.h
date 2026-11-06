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
    /// Class for EStep (RS)
    class EStep {
        ReconstructionModel *reconstructor;
        Array<double>& slice_potential;

    public:
        EStep(ReconstructionModel *reconstructor,
            Array<double>& slice_potential) :
            reconstructor(reconstructor),
            slice_potential(slice_potential) {}

        void operator()(const blocked_range<size_t>& r) const {
            RealImage slice;

            for (size_t inputIndex = r.begin(); inputIndex < r.end(); inputIndex++) {
                // read the current slice

                if (reconstructor->_no_masking_background)
                    slice = reconstructor->_not_masked_slices[inputIndex];
                else
                    slice = reconstructor->_slices[inputIndex];

                //read current weight image
                RealImage& weight = reconstructor->_weights[inputIndex];
                memset(weight.Data(), 0, sizeof(RealPixel) * weight.NumberOfVoxels());

                double num = 0;
                //Calculate error, voxel weights, and slice potential
                for (int l = 0; l < slice.Attributes().T(); l++ )
                    for (size_t i = 0; i < reconstructor->_volcoeffs[inputIndex].size(); i++)
                        for (size_t j = 0; j < reconstructor->_volcoeffs[inputIndex][i].size(); j++)
                            num += calculateSlice(slice, weight, inputIndex, i, j, l);

                //evaluate slice potential
                if (num > 0)
                    slice_potential[inputIndex] = sqrt(slice_potential[inputIndex] / num);
                else
                    slice_potential[inputIndex] = -1; // slice has no unpadded voxels
            }
        }

        int calculateSlice(  RealImage& slice,
                                    RealImage& weight,
                                    int inputIndex, 
                                    int i, 
                                    int j, 
                                    int t ) const {
            if (reconstructor->InSliceMask(inputIndex,i, j)) {
                //bias correct and scale the slice
                slice(i, j, 0, t) *= exp(-reconstructor->_bias[inputIndex](i, j, 0)) * reconstructor->_scale[inputIndex];

                //number of volumetric voxels to which
                // current slice voxel contributes
                const size_t n = reconstructor->_volcoeffs[inputIndex][i][j].size();

                // if n == 0, slice voxel has no overlap with volumetric ROI, do not process it

                if (n > 0 && reconstructor->_simulated_weights[inputIndex](i, j, 0, t) > 0) {
                    slice(i, j, 0, t) -= reconstructor->_simulated_slices[inputIndex](i, j, 0, t);

                    //calculate norm and voxel-wise weights

                    //Gaussian distribution for inliers (likelihood)
                    const double g = reconstructor->G(slice(i, j, 0, t), reconstructor->_sigma);
                    //Uniform distribution for outliers (likelihood)
                    const double m = reconstructor->M(reconstructor->_m);

                    //voxel_wise posterior
                    const double w = g * reconstructor->_mix / (g * reconstructor->_mix + m * (1 - reconstructor->_mix));
                    weight.PutAsDouble(i, j, 0, t, w);

                    //calculate slice potentials
                    if (reconstructor->_simulated_weights[inputIndex](i, j, 0, t) > 0.99) {
                        slice_potential[inputIndex] += (1 - w) * (1 - w);
                        return 1;
                    }
                } else
                    weight.PutAsDouble(i, j, 0, t, 0);
            }
            return 0;
        }

        void operator()() const {
            parallel_for(blocked_range<size_t>(0, reconstructor->_slices.size()), *this);
        }
    };

}