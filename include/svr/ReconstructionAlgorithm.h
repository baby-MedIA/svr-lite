// Copyright © King's College London
// SPDX-License-Identifier: Apache-2.0
//
// This file has been modified from the original source.

#pragma once

// SVR
#include "svr/Application.h"
#include "svr/ReconstructionModel.h"
#include "svr/ReconstructionParameters.h"
#include <filesystem>

using namespace mirtk;

namespace svr {
    /***
     * Reconstruction algorithm
     * Defines the logical steps of the reconstruction process.
     */
    class ReconstructionAlgorithm : public Application {
        public:
            ReconstructionAlgorithm(ReconstructionParameters &inputData);

            // Run reconstruction algorithm
            int Run();

        private:

            ReconstructionModel _reconstruction;
            ReconstructionParameters& _parameters;

            Array<RealImage> _stacks;
            Array<RigidTransformation> _stackTransformations; // recon _transformations (check that these are per stack though)
            RigidTransformation _transformationReconToRef;

            Array<double> _meanDisplacement; // Mean Displacement
            Array<double> _meanWeightedDisplacement;
            Array<Array<double>> _entropy;

            /*************************
             * Main algorithm stages *
             *************************/

            /***
             * Initialise data structures ahead of reconstruction
             */
            int Initialise();

            /***
             * Pre-processing steps ahead of recontructions
             */
            int PreProcess();

            /***
             * Process the 4D reconstruction algoritm
             */
            int Process();

            /*** 
             * Evaluate
             */
            void Evaluate(const boost::timer::cpu_timer& reconstruction_timer);

            /***
             * Post processing, save outputs and clean up
             */
            int PostProcess();


            /************************
             * Supporting Functions *
             ************************/

            void CalculateSliceThickness();

            void ReadPackagesParam();

            void ConfigureDebugMode();

            void ReadStackFiles();

            void ReadStackTransformationFiles();

            void CropStacks();

            void RemoveBlackBackground();

            void InitialiseVolume(mirtk::RealImage &maskCropped);

            void SetupVolumeMask(std::unique_ptr<mirtk::RealImage> &reconstructionMask, mirtk::RealImage &maskCropped);

            void SplitStacksIntoPackages();

            void CalculateDisplacements(int &iter, RealImage& reconstructed);

            void SaveIntermediateSimulatedSlices(int iter, int i);

            int CalcReconstuctionIterations(int iter);

            void SaveIntermediateSimulatedSlices(int iter);

            void SetSmoothingParameters(int iter);

            void SliceToVolumeRegistration(int iter);


    };
}