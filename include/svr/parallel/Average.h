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

    /// Class for computing average of all stacks
    class Average {
        const ReconstructionModel *reconstructor;
        const Array<RealImage>& stacks;
        const Array<RigidTransformation>& stack_transformations;

        /// Padding value in target (voxels in the target image with this value will be ignored)
        double targetPadding;

        /// Padding value in source (voxels outside the source image will be set to this value)
        double sourcePadding;

        double background;

        // Volumetric registrations are stack-to-template while slice-to-volume
        // registrations are actually performed as volume-to-slice
        // (reasons: technicalities of implementation) so transformations need to be inverted beforehand.

        bool linear;

    public:
        RealImage average;
        RealImage weights;

        Average(
            const ReconstructionModel *reconstructor,
            const Array<RealImage>& stacks,
            const Array<RigidTransformation>& stack_transformations,
            double targetPadding,
            double sourcePadding,
            double background,
            bool linear = false) :
            reconstructor(reconstructor),
            stacks(stacks),
            stack_transformations(stack_transformations),
            targetPadding(targetPadding),
            sourcePadding(sourcePadding),
            background(background),
            linear(linear) {
            average.Initialize(reconstructor->_reconstructed.Attributes());
            weights.Initialize(reconstructor->_reconstructed.Attributes());
        }

        Average(Average& x, split) : Average(x.reconstructor, x.stacks, x.stack_transformations, x.targetPadding, x.sourcePadding, x.background, x.linear) {}

        void operator()(const blocked_range<size_t>& r) {
            GenericLinearInterpolateImageFunction<RealImage> interpolator;
            ImageTransformation imagetransformation;
            RealImage image;

            for (size_t i = r.begin(); i < r.end(); i++) {
                image.Initialize(reconstructor->_reconstructed.Attributes());

                imagetransformation.Input(&stacks[i]);
                imagetransformation.Transformation(&stack_transformations[i]);
                imagetransformation.Output(&image);
                imagetransformation.TargetPaddingValue(targetPadding);
                imagetransformation.SourcePaddingValue(sourcePadding);
                imagetransformation.Interpolator(&interpolator);
                imagetransformation.Run();

                const RealPixel *pi = image.Data();
                RealPixel *pa = average.Data();
                RealPixel *pw = weights.Data();
                for (int p = 0; p < average.NumberOfVoxels(); p++) {
                    if (pi[p] != background) {
                        pa[p] += pi[p];
                        pw[p]++;
                    }
                }
            }
        }

        void join(const Average& y) {
            average += y.average;
            weights += y.weights;
        }

        void operator()() {
            parallel_reduce(blocked_range<size_t>(0, stacks.size()), *this);
        }
    };
}