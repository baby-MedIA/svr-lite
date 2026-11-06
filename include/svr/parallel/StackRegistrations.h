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

    // class for global 3D stack registration
    class StackRegistrations {
        const ReconstructionModel *reconstructor;
        const Array<RealImage>& stacks;
        Array<RigidTransformation>& stack_transformations;
        const RealImage& target;
        const RigidTransformation& offset;

    public:
        StackRegistrations(
            const ReconstructionModel *reconstructor,
            const Array<RealImage>& stacks,
            Array<RigidTransformation>& stack_transformations,
            const RealImage& target,
            const RigidTransformation& offset) :
            reconstructor(reconstructor),
            stacks(stacks),
            stack_transformations(stack_transformations),
            target(target),
            offset(offset) {}

        void operator()(const blocked_range<size_t>& r) const {
            RealImage r_source, r_target;
            GenericLinearInterpolateImageFunction<RealImage> interpolator;
            ResamplingWithPadding<RealPixel> resampling(1.2, 1.2, 1.2, 0);
            resampling.Interpolator(&interpolator);

            ParameterList params;
            Insert(params, "Transformation model", "Rigid");
            if (reconstructor->_nmi_bins > 0)
                Insert(params, "No. of bins", reconstructor->_nmi_bins);

            if (reconstructor->_masked_stacks)
                Insert(params, "Background value", 0);
            else
                Insert(params, "Background value for image 1", 0);

            // Insert(params, "Background value", 0);
            // Insert(params, "Image (dis-)similarity measure", "NCC");
            // Insert(params, "Image interpolation mode", "Linear");
            // string type = "box";
            // string units = "mm";
            // double width = 0;
            // Insert(params, string("Local window size [") + type + string("]"), ToString(width) + units);

            Transformation *dofout;
            GenericRegistrationFilter registration;
            registration.Parameter(params);
            registration.Output(&dofout);

            for (size_t i = r.begin(); i != r.end(); i++) {

                r_source.Initialize(stacks[i].Attributes());
                resampling.Input(&stacks[i]);
                resampling.Output(&r_source);
                resampling.Run();

                r_target.Initialize(target.Attributes());
                resampling.Input(&target);
                resampling.Output(&r_target);
                resampling.Run();

                const Matrix& mo = offset.GetMatrix();
                stack_transformations[i].PutMatrix(stack_transformations[i].GetMatrix() * mo);

                registration.Input(&r_target, &r_source);
                registration.InitialGuess(&stack_transformations[i]);
                registration.GuessParameter();
                registration.Run();

                unique_ptr<RigidTransformation> rigidTransf(dynamic_cast<RigidTransformation*>(dofout));
                stack_transformations[i] = *rigidTransf;

                stack_transformations[i].PutMatrix(stack_transformations[i].GetMatrix() * mo.Inverse());
            }
        }

        void operator()() const {
            parallel_for(blocked_range<size_t>(0, stacks.size()), *this);
        }
    };
}