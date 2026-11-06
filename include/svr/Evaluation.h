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
     * Reconstruction evaluation tool
     * Defines the logical steps of the evaluation process.
     */
    class Evaluation : public Application {
        public:
            // Constructor
            Evaluation(const bool debug, const bool verbose, const bool profile);

            // Run reconstruction algorithm
            int run(const string& stackFile,
                    double& thickness,
                    const string& maskFile,
                    const string& reconVolume,
                    const bool stackIntensityMatching,
                    const int& iterations);


            // Read stacks from file into stack array
            void ReadStack(const string& stackFile, RealImage& stack, int& tSize);

            // Read masks from file into mask array
            void ReadMask(const string& maskFile, RealImage& mask);

            // Calculates thickness of slices in stacks respectively
            void CalculateSliceThickness(RealImage &stacks, double &thickness);

            // Read volume from file
            void ReadVolume(const string &volumeFile, RealImage &volume, int &tSize);

            // Create mask from volume padding
            void ExtractVolumeMask(const RealImage& volume, RealImage& mask);

            // Similarity Metrics:

            // Calcuale Normalised Cross-Correlation
            double CalculateNCC(const RealImage& target, 
                                       const RealImage& source, 
                                       const RealImage& masks,
                                      const RealImage& bias,
                                      const vector<double>& scale);

            // Calcuale Normalised Root Mean Squared Error
            double CalculateNRMSE(const RealImage& target, 
                                        const RealImage& source, 
                                        const RealImage& masks,
                                      const RealImage& bias,
                                      const vector<double>& scale);
    };
}