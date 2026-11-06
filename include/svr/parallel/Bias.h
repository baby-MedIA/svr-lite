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
    class Bias {
        ReconstructionModel *reconstructor;

    public:
        Bias(ReconstructionModel *reconstructor) : reconstructor(reconstructor) {}

        void operator()(const blocked_range<size_t>& r) const {
            RealImage slice, wb, wresidual;
            GaussianBlurring<RealPixel> gb(reconstructor->_sigma_bias);

            for (size_t inputIndex = r.begin(); inputIndex < r.end(); inputIndex++) {
                // read the current slice
                slice = reconstructor->_slices[inputIndex];

                //alias the current bias image
                RealImage& b = reconstructor->_bias[inputIndex];

                //prepare weight image for bias field
                wb = reconstructor->_weights[inputIndex];

                wresidual.Initialize(slice.Attributes());

                for (int l = 0; l < slice.GetT(); l++)
                    for (int j = 0; j < slice.GetY(); j++)
                        for (int i = 0; i < slice.GetX(); i++)
                            if (reconstructor->InSliceMask(inputIndex,i, j))
                            {
                                if (reconstructor->_simulated_weights[inputIndex](i, j, 0, l) > 0.99)
                                {
                                    // bias-correct and scale current slice
                                    const double eb = exp(-b(i, j, 0));
                                    slice(i, j, 0, l) *= eb * reconstructor->_scale[inputIndex];

                                    // calculate weight image
                                    wb(i, j, 0, l) *= slice(i, j, 0, l);

                                    // calculate weighted residual image make sure it is far from zero 
                                    // to avoid numerical instability
                                    if ((reconstructor->_simulated_slices[inputIndex](i, j, 0, l) > 1) 
                                        && (slice(i, j, 0, l) > 1))
                                        wresidual(i, j, 0, l) = log(slice(i, j, 0, l) / 
                                            reconstructor->_simulated_slices[inputIndex](i, j, 0, l)) * wb(i, j, 0, l);
                                }
                                else
                                {
                                    // do not take into account this voxel when calculating bias field
                                    wresidual(i, j, 0, l) = 0;
                                    wb(i, j, 0, l) = 0;
                                }
                            }

                //calculate bias field for this slice
                //smooth weighted residual
                gb.Input(&wresidual);
                gb.Output(&wresidual);
                gb.Run();

                //smooth weight image
                gb.Input(&wb);
                gb.Output(&wb);
                gb.Run();

                //update bias field
                double sum = 0;
                double num = 0;

                RealImage b_init = b; // make a copy of b before incrementing time values.

                for (int l = 0; l < slice.GetT(); l++)
                    for (int i = 0; i < slice.GetX(); i++)
                        for (int j = 0; j < slice.GetY(); j++)
                            if (reconstructor->InSliceMask(inputIndex,i, j))
                            {
                                float biasVal = 0;
                                if (wb(i, j, 0, l) > 0) {
                                    biasVal = wresidual(i, j, 0, l) / wb(i, j, 0, l) / slice.GetT(); // add mean bias across T.
                                    b(i, j, 0) += biasVal;
                                }

                                // Must use original value b_init(i,j) of b(i,j) when calculating sum for mean 
                                // so that biasVal at previous iterations of t isn't counted cumulatively.
                                // (a side effect of reducing from 4D to 3D within this calculation).
                                sum += b_init(i, j, 0) + biasVal; 

                                num++;
                            }

                //normalize bias field to have zero mean
                if (!reconstructor->_global_bias_correction && num > 0) {
                    const double mean = sum / num;
                    for (int i = 0; i < slice.GetX(); i++)
                        for (int j = 0; j < slice.GetY(); j++)
                            if (reconstructor->InSliceMask(inputIndex,i, j))
                                b(i, j, 0) -= mean;
                }
            }
        }

        void operator()() const {
            parallel_for(blocked_range<size_t>(0, reconstructor->_slices.size()), *this);
        }
    };

}