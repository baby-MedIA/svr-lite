// Copyright © King's College London
// SPDX-License-Identifier: Apache-2.0
//
// This file has been modified from the original source.

#pragma once

// SVR
#include "svr/Common.h"

namespace svr {

    /**
     * Intensity Matching mode options
     * 
     * NO_RESCALING = Do not apply intensity matching.
     * RESCALE_STACKS_ONLY = Apply intensity matching to stacks only.
     * FULL_RESCALING = Apply intensity matching to stacks and slices throughout.
     */
    enum IM_MODE { NO_RESCALING=0, RESCALE_STACKS_ONLY, FULL_RESCALE_EXCEPT_STACKS, FULL_RESCALING};

    /**
     * Slice masking mode options
     * 
     * NO_SLICE_MASKS = Do not create slice masks or appy them anywhere. (Original 4D algorithm design)
     * MASK_SLICES_ONLY_SVR = Use slice masks during slice-to-volume registration only.
     * MASK_SLICES_NOT_SVR = Use slice masks everywhere except slice-to-volume registration. (Original 3D algorithm design - see legacy note below)
     * MASK_SLICES_ALWAYS = Use slice masks everywhere.
     * LEGACY_APPLY_SLICE_MASKS = TEMPORARY. Legacy 3D behaviour for regression testing. Apply slice masks to slices. TO BE RETIRED in favour of MASK_SLICES_NOT_SVR.
     */
    enum MASK_MODE { NO_SLICE_MASKS=0, MASK_SLICES_ONLY_SVR, MASK_SLICES_NOT_SVR, MASK_SLICES_ALWAYS };

    /***
     * Input parameters for Reconstruction 4D
     */
    struct ReconstructionParameters { 
        string outputName;                  // Output file name

        //Files
        Array<string> stackFiles;           // The input stacks (Nifti format)
        vector<string> stackRegFolder;          // The transformations of the input stack to template in \'dof\' format used in IRTK. Only rough alignment with correct orientation and some overlap is needed. Use \'id\' for an identity transformation for at leastone stack. The first stack with \'id\' transformation will be resampled as template.
        string frameRegfolder;                      // Folder for slice-to-volume registrations, if given
        string sliceLocRegFolder;  // Folder for slice-location registrations, if given
        string maskFile;                    // Paths of mask file (single)
        string logID = "";                  // Prefix for the log file.
        string infoFilename = "info.tsv";   // File name for slice information in tab-separated columns.

        // Data processing parameters
        int iterations = 4;                 // Number of registration iterations
        int srIterations = 21;              // Number of SR reconstruction iterations for final iteration
        int templateNumber = 0;             // Stack number of target for stack-to-stack registration.
        Array<double> thickness;            // Give slice thickness. 
        int nmiBins = 16;                   // Numbers of NMI bins for registration
        Array<int> forceExcludedSlices;     // Force exclusion of image-slices with these indices.
        Array<int> forceExcludedStacks;     // Force exclusion of stacks
        vector<int> packages;               // Array of number of packages for each stack
        double resolution = 0.75;           // Isotropic resolution of the volume.
        double timeResolution = 0;          // Resolution on t-axis in milliseconds. If 0 use average of stacks.
        double sigma = 20;                  // Stdev for bias field.
        double averageValue = 700;          // Average intensity value for stacks.
        double delta = 150;                 // Parameter to define what is an edge. 
        double lambda = 0.02;               // Smoothing parameter.
        double lastIterLambda = 0.01;       // Smoothing parameter for last iteration. 
        int levels = 3;                     // Multiresolution smoothing with given number of levels.
        double smoothMask = 4;              // Smooth the mask to reduce artefacts of manual segmentation.

        int intensityMatchingMode = IM_MODE::FULL_RESCALING;
        int sliceMaskingMode = MASK_MODE::NO_SLICE_MASKS;

        // Flags
        bool removeBlackBackground = false; // Generate mask from black background
        bool structural = false;            // Flag for struture-based exclusion of slices
        bool svrMaskSlice = false;          // Mask slices for slice-to-volume registration step.
        bool robustStatistics = true;       // Flag to swich the robust statistics on and off
        bool robustSlicesOnly = false;      // Robust statistics for exclusion of slices only
        bool saveSlicesFlag = false;        // Flag for saving slices
        bool saveEvaluationData = false;    // Flag for saving evaluation data
        bool debug = false;                 // Debug mode
        bool verbose = false;               // Verbose debugging mode
        bool profile = false;               // Profiling mode

        // LEGACY OPTIONS FOR REGRESSION TESTING ONLY
        bool noStackZCrop = false;          // Crop stacks on the Z axis.
        bool performSecondStackReg = false; // Perform second stack registration
        bool lastIterRSCalcs = false;       // Calculate robust statistics on last iteration

        public:

        // Validate set parameters
        bool Validate();
    };

    ostream& operator<<(ostream& os, const ReconstructionParameters& params);  
}