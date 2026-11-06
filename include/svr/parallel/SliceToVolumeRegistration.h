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

    /***
     * Slice to Volume registration
     * Registers each individual slice to the reconstructed volume.
     * 
     * Primary inputs for registration:
     * Target = Original slice at first timepoint.
     * Source = 3D volume extracted from first timepoint.
     * Initial Guess = Stack-to-stack transformation or the transformation from a previous iteration.
     * Transformation model = Rigid
     * Image (dis-)similarity measure = NMI
     * Background value for image 2 = -1
     * 
     * If slice masking is set then slice mask is written to local slice 
     * copy before registration and target background set to 0.
     */
    class SliceToVolumeRegistration {
        ReconstructionModel *reconstructor;
        int iteration = 0;
        bool maskSlices = false;

        const int SOURCE_PADDING = -1;


    public:
        SliceToVolumeRegistration(ReconstructionModel *reconstructor, 
                                    int iter, bool maskSlices) 
            : reconstructor(reconstructor),
              iteration(iter),
              maskSlices(maskSlices)
              {}

        void operator()(const blocked_range<size_t>& r) const {
                
            const ImageAttributes& attr = reconstructor->_reconstructed.Attributes();

            // Convert source to GreyImage for efficient registration.
            GreyImage source = reconstructor->_reconstructed.GetRegion(0, 0, 0, 0, attr._x, attr._y, attr._z, 1);   

            ParameterList params;
            Insert(params, "Transformation model", "Rigid");
            Insert(params, "Image (dis-)similarity measure", "NMI");
            if (reconstructor->_nmi_bins > 0)
                Insert(params, "No. of bins", reconstructor->_nmi_bins);

            Insert(params, "Background value for image 2", SOURCE_PADDING);

            Transformation *dofout;
            GenericRegistrationFilter registration;
            registration.Parameter(params);
            registration.Output(&dofout);

            for (size_t inputIndex = r.begin(); inputIndex != r.end(); inputIndex++) {
                if (reconstructor->_slice_excluded[inputIndex])
                    continue;

                // get pixel value min and max
                const ImageAttributes& target_attr = reconstructor->_slices[inputIndex].Attributes(); 

                RealImage targetReal = reconstructor->_slices[inputIndex]
                                        .GetRegion(0,0,0,0,target_attr._x, target_attr._y, target_attr._z, 1);

                if (maskSlices) {
                    const RealImage& mask = reconstructor->_slice_masks[inputIndex];
                    targetReal *= mask;
                    Insert(params, "Background value for image 1", 0);
                }

                // Convert target to GreyImage for efficient registration.
                GreyImage target = targetReal;
                
                // SOURCE
                GreyPixel smin, smax;
                target.GetMinMax(&smin, &smax);
                if (smax > 0 && (smax - smin) > 1) {
                    // put origin to zero
                    RigidTransformation offset;
                    ResetOrigin(target, offset);

                   const Matrix& mo = offset.GetMatrix();
                    auto& transformation = reconstructor->_transformations[inputIndex];
                    transformation.PutMatrix(transformation.GetMatrix() * mo);
                    registration.Input(&target, &source);
                    registration.InitialGuess(&transformation);
                    registration.GuessParameter();
                    registration.Run();
                    unique_ptr<RigidTransformation> rigidTransf(dynamic_cast<RigidTransformation*>(dofout));
                    transformation = *rigidTransf;

                    //undo the offset
                    transformation.PutMatrix(transformation.GetMatrix() * mo.Inverse());
                }
            }
        }

        void operator()() const {
            logger::trace("Running slice to volume registration");
            parallel_for(blocked_range<size_t>(0, reconstructor->_slices.size()), *this);
            logger::trace("Completed slice to volume registration");
        }
    };
}