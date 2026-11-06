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

    /// Class for MStep (RS)
    class MStep {
        ReconstructionModel *reconstructor;

    public:
        double sigma;
        double mix;
        double num;
        double min;
        double max;

        MStep(ReconstructionModel *reconstructor) : reconstructor(reconstructor) {
            sigma = 0;
            mix = 0;
            num = 0;
            min = voxel_limits<RealPixel>::max();
            max = voxel_limits<RealPixel>::min();
        }

        MStep(MStep& x, split) : MStep(x.reconstructor) {}

        void operator()(const blocked_range<size_t>& r) {
            RealImage slice;

            for (size_t inputIndex = r.begin(); inputIndex < r.end(); inputIndex++) {
                // read the current slice
                slice = reconstructor->_slices[inputIndex];

                //calculate error
                for (int l = 0; l < slice.GetT(); l++)
                    for (int j = 0; j < slice.GetY(); j++)
                        for (int i = 0; i < slice.GetX(); i++)
                            CalculateError(i, j, l, inputIndex, slice);
            } //end of loop for a slice inputIndex
        }

        void CalculateError(int &i, int &j, int l, size_t &inputIndex, RealImage &slice)
        {
            if (reconstructor->InSliceMask(inputIndex,i, j)) 
            {
                // bias correct and scale the slice
                slice(i, j, 0, l) *= exp(-reconstructor->_bias[inputIndex](i, j, 0)) 
                                        * reconstructor->_scale[inputIndex];

                // otherwise the error has no meaning - it is equal to slice intensity
                if (reconstructor->_simulated_weights[inputIndex](i, j, 0, l) > 0.99)
                {
                    slice(i, j, 0, l) -= reconstructor->_simulated_slices[inputIndex](i, j, 0, l);

                    // sigma and mix
                    const double e = slice(i, j, 0, l);
                    sigma += e * e * reconstructor->_weights[inputIndex](i, j, 0, l);
                    mix += reconstructor->_weights[inputIndex](i, j, 0, l);

                    //_m
                    if (e < min)
                        min = e;
                    if (e > max)
                        max = e;

                    num++;
                }
            }
        }


        void join(const MStep& y) {
            if (y.min < min)
                min = y.min;
            if (y.max > max)
                max = y.max;

            sigma += y.sigma;
            mix += y.mix;
            num += y.num;
        }

        void operator()() {
            parallel_reduce(blocked_range<size_t>(0, reconstructor->_slices.size()), *this);
        }
    };
}