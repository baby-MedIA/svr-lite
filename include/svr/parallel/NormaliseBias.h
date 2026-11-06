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

    class NormaliseBias {
        ReconstructionModel *reconstructor;

    public:
        RealImage bias;
        RealImage volweight3d;

        NormaliseBias(ReconstructionModel *reconstructor) : reconstructor(reconstructor) {
            ImageAttributes attr = reconstructor->_reconstructed.Attributes();
            attr._t = 1;
            bias.Initialize(attr);
            volweight3d.Initialize(attr);
        }

        NormaliseBias(NormaliseBias& x, split) : NormaliseBias(x.reconstructor) {}

        void operator()(const blocked_range<size_t>& r) {
            RealImage b;
            for (size_t inputIndex = r.begin(); inputIndex < r.end(); inputIndex++) {
                if (reconstructor->_volcoeffs[inputIndex].empty())
                    continue;

                // alias to the current slice
                const RealImage& slice = reconstructor->_slices[inputIndex];

                //read the current bias image
                b = reconstructor->_bias[inputIndex];

                //read current scale factor
                const double scale = reconstructor->_scale[inputIndex];


                    // Only execute on t=0 for bias.
                    for (int j = 0; j < slice.GetY(); ++j)
                        for (int i = 0; i < slice.GetX(); ++i) 
                            if (reconstructor->InSliceMask(inputIndex,i, j) && scale > 0)
                                b(i, j, 0) -= log(scale);   

                //Distribute slice intensities to the volume
                for (int l = 0; l < slice.GetT(); l++)
                    for (size_t i = 0; i < reconstructor->_volcoeffs[inputIndex].size(); i++)
                        for (size_t j = 0; j < reconstructor->_volcoeffs[inputIndex][i].size(); j++)
                            if (reconstructor->InSliceMask(inputIndex,i, j)) {
                                //add contribution of current slice voxel to all voxel volumes
                                //to which it contributes
                                for (size_t k = 0; k < reconstructor->_volcoeffs[inputIndex][i][j].size(); k++) {
                                    const POINT3D& p = reconstructor->_volcoeffs[inputIndex][i][j][k];

                                    if (l==0) //only calculate bias at t=0
                                        bias(p.x, p.y, p.z) += p.value * b(i, j, 0);
                                    volweight3d(p.x, p.y, p.z) += p.value;
                                }
                            }
                //end of loop for a slice inputIndex
            }
        }

        void join(const NormaliseBias& y) {
            bias += y.bias;
            volweight3d += y.volweight3d;
        }

        void operator()() {
            parallel_reduce(blocked_range<size_t>(0, reconstructor->_slices.size()), *this);
        }
    };

}
