// Copyright © King's College London
// SPDX-License-Identifier: Apache-2.0
//
// This file has been modified from the original source.


// SVR
#include "svr/ReconstructionParameters.h"
#include <ostream>
#include <filesystem>

using namespace std;


namespace svr {

    ostream& operator<<(ostream& os, const ReconstructionParameters& params) {
        os << "outputName: " << params.outputName << endl;
        os << "stackFiles: "; for (auto n: params.stackFiles) { os << n << ", "; } os << endl;
        os << "templateNumber: " << params.templateNumber << " (Not zero indexed=" << (params.templateNumber+1) << ")" << endl;
        os << "packages: "; for (auto n: params.packages) { os << n << ", "; } os << endl;
        os << "dofinPaths: "; for (auto n: params.stackRegFolder) { os << n << ", "; } os << endl;
        os << "thickness: "; for (auto n: params.thickness) { os << n << ", "; } os << endl;
        os << "folder: " << params.frameRegfolder << endl;
        os << "sliceTransformationsFolder: " << params.sliceLocRegFolder << endl;
        os << "resolution: " << params.resolution << endl;
        os << "time_resolution: " << params.timeResolution << endl;
        os << "iterations: " << params.iterations << endl;
        os << "srIterations: " << params.srIterations << endl;
        os << "sigma: " << params.sigma << endl;
        os << "averageValue: " << params.averageValue << endl;
        os << "delta: " << params.delta << endl;
        os << "lambda: " << params.lambda << endl;
        os << "lastIterLambda: " << params.lastIterLambda << endl;
        os << "levels: " << params.levels << endl;
        os << "smoothMask: " << params.smoothMask << endl;
        os << "forceExcludedSlices: "; for (auto n: params.forceExcludedSlices) { os << n << ", "; } os << endl;
        os << "robustStatistics: " << params.robustStatistics << endl;
        os << "robustSlicesOnly: " << params.robustSlicesOnly << endl;
        os << "structural: " << params.structural << endl;
        os << "maskFile: " << params.maskFile << endl;
        os << "removeBlackBackground: " << params.removeBlackBackground << endl;
        os << "nmiBins: " << params.nmiBins << endl;
        os << "logID: " << params.logID << endl;
        os << "infoFilename: " << params.infoFilename << endl;

        os << "intensityMatchingMode: " << params.intensityMatchingMode << endl;
        os << "sliceMaskingMode: " << params.sliceMaskingMode << endl;
        os << "LEGACY: noStackZCrop: " << params.noStackZCrop << endl;
        os << "LEGACY: performSecondStackReg: " << params.performSecondStackReg << endl;
        os << "LEGACY: lastIterRSCalcs: " << params.lastIterRSCalcs << endl;

        os << "saveSlicesFlag: " << params.saveSlicesFlag << endl;
        os << "debug: " << params.debug << endl;
        os << "verbose: " << params.verbose << endl;
        os << "profile: " << params.profile << endl;
        return os;
    }

    bool ReconstructionParameters::Validate() {
        if (!stackRegFolder.empty() && stackRegFolder.size() < stackFiles.size()) {
            logger::error("Count of dof files should equal to stack count!");
            return false;
        }

        if (maskFile.empty() && !removeBlackBackground) {
            logger::error("No mask file and removeBlackBackground not specified. A mask is required.");
            return false;
        }

        // check that conflicting transformation folders haven't been given
        if (!frameRegfolder.empty() && !sliceLocRegFolder.empty()) {
            logger::error("Can not use both -transformations and -slice_transformations arguments.");
            return 1;
        }

        // output name
        if(outputName.empty()){
            logger::error("Output name must be specified");
            return 1;
        }

        // stack files
        if(stackFiles.size() < 1){
            logger::error("Must provide at least 1 stack file");
            return 1;
        }

        for(string fileName : stackFiles)
            if (!std::filesystem::exists(fileName)){
                logger::error("Invalid stack filename: " + fileName);
                return 1;
            }

        // dofin paths
        if (!stackRegFolder.empty() && stackRegFolder.size() < stackFiles.size()){
            logger::error("Count of masks should equal to stack count!");
            return 1;
        }

        for(string fileName : stackRegFolder)
            if (!std::filesystem::exists(fileName)){
                logger::error("Invalid dofin filename: " + fileName);
                return 1;
            }

        // thickness
        if (!thickness.empty() && thickness.size() < stackFiles.size()){
            logger::error("Count of thickness values should equal to stack count!");
            return 1;
        }

        // maskFile
        if (!maskFile.empty() && !std::filesystem::exists(maskFile)){
            logger::error("Invalid mask filename: " + maskFile);
            return 1;
        }
        
        // Check whether the template stack can be identified
        if (templateNumber < 0) {
            logger::error("Please identify the template by assigning id transformation.");
            return 1;
        }

        return true;
    }
}