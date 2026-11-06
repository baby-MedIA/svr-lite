// Copyright © King's College London
// SPDX-License-Identifier: Apache-2.0
//
// This file has been modified from the original source.

#pragma once

// SVR
#include "svr/ReconstructionModel.h"
#include "mirtk/VoxelFunction.h"

using namespace mirtk;
using namespace svr;
using namespace svr::Utility;

namespace svr::Parallel {

    class AdaptiveRegularizationStepA : public VoxelFunction{
        const ReconstructionModel *reconstructor;
        Array<RealImage>& b;
        const Array<double>& factor;
        const RealImage& original;

    public:
        AdaptiveRegularizationStepA(
            const ReconstructionModel *reconstructor,
            Array<RealImage>& b,
            const Array<double>& factor,
            const RealImage& original) :
            reconstructor(reconstructor),
            b(b),
            factor(factor),
            original(original) {}

        void operator()(int x, int y, int z, int t, realt *p) const {
            for (size_t i = 0; i < 13; i++) {
                const int xx = x + reconstructor->_directions[i][0];
                const int yy = y + reconstructor->_directions[i][1];
                const int zz = z + reconstructor->_directions[i][2];
                if ((xx >= 0) && (xx < reconstructor->_reconstructed.GetX()) && 
                    (yy >= 0) && (yy < reconstructor->_reconstructed.GetY()) && 
                    (zz >= 0) && (zz < reconstructor->_reconstructed.GetZ()) && 
                    (reconstructor->_confidence_map(x, y, z, t) > 0) && 
                    (reconstructor->_confidence_map(xx, yy, zz, t) > 0)) 
                    // (reconstructor->_mask(x, y, z) > 0) && 
                    // (reconstructor->_mask(xx, yy, zz) > 0))
                {  
                    const double diff = (original(xx, yy, zz, t) - original(x, y, z, t)) * sqrt(factor[i]) / reconstructor->_delta;
                    b[i](x, y, z, t) = factor[i] / sqrt(1 + diff * diff);
                } else
                    b[i](x, y, z, t) = 0;
            }
        }
    };

    class AdaptiveRegularizationStepB : public VoxelFunction{
        ReconstructionModel *reconstructor;
        const Array<RealImage>& b;
        const RealImage& original;

    public:
        AdaptiveRegularizationStepB(
            ReconstructionModel *reconstructor,
            Array<RealImage>& b,
            const RealImage& original) :
            reconstructor(reconstructor),
            b(b),
            original(original) {}

        void operator()(int x, int y, int z, int t, realt *p) const { 
            if (reconstructor->_confidence_map(x, y, z, t) > 0) { 
            //if (reconstructor->_mask(x, y, z) > 0) { //Future potential improvement
                double val = 0;
                double sum = 0;
                for (int i = 0; i < 13; i++) {
                    const int xx = x + reconstructor->_directions[i][0];
                    const int yy = y + reconstructor->_directions[i][1];
                    const int zz = z + reconstructor->_directions[i][2];
                    if ((xx >= 0) && (xx < reconstructor->_reconstructed.GetX()) && 
                        (yy >= 0) && (yy < reconstructor->_reconstructed.GetY()) && 
                        (zz >= 0) && (zz < reconstructor->_reconstructed.GetZ()) && 
                        (reconstructor->_confidence_map(xx, yy, zz, t) > 0)) {
                        // (reconstructor->_mask(xx, yy, zz) > 0)) {//Future potential improvement
                        val += b[i](x, y, z, t) * original(xx, yy, zz, t);
                        sum += b[i](x, y, z, t);
                    }
                }

                for (int i = 0; i < 13; i++) {
                    const int xx = x - reconstructor->_directions[i][0];
                    const int yy = y - reconstructor->_directions[i][1];
                    const int zz = z - reconstructor->_directions[i][2];
                    if ((xx >= 0) && (xx < reconstructor->_reconstructed.GetX()) && 
                        (yy >= 0) && (yy < reconstructor->_reconstructed.GetY()) && 
                        (zz >= 0) && (zz < reconstructor->_reconstructed.GetZ()) && 
                        (reconstructor->_confidence_map(xx, yy, zz, t) > 0)) {
                        //(reconstructor->_mask(xx, yy, zz) > 0)) {//Future potential improvement
                            val += b[i](x, y, z, t) * original(xx, yy, zz, t);
                            sum += b[i](x, y, z, t);
                        }
                }

                val -= sum * original(x, y, z, t);
                val = original(x, y, z, t) + reconstructor->_alpha * reconstructor->_lambda / (reconstructor->_delta * reconstructor->_delta) * val;
                reconstructor->_reconstructed(x, y, z, t) = val;
            }
        }
    };
}