// Copyright © King's College London
// SPDX-License-Identifier: Apache-2.0
//
// This file has been modified from the original source.

// SVR
#include "svr/ReconstructionAlgorithm.h"


#define SVR_TOOL
#include "svr/Profiling.h"

using namespace std;
using namespace mirtk;
using namespace svr;
using namespace boost::program_options;

namespace reconstruct_tool {

    options_description getOptions(ReconstructionParameters& parameters) {
        options_description opts("Options");
        opts.add_options()
            ("mask", value<string>(&parameters.maskFile), "Binary mask to define the region of interest. [Required]")
            ("packages", value<vector<int>>(&parameters.packages)->multitoken(), "Number of packages used during acquisition per stack, or default number for all stacks. Stacks will be split into packages.")
            ("target_stack", value<int>(&parameters.templateNumber), "Stack number of target for stack-to-stack registration.")
            ("dofin", value<vector<string>>(&parameters.stackRegFolder)->multitoken(), "The transformations of the input stack to template in \'dof\' format used in IRTK. Only rough alignment with correct orientation and some overlap is needed. Use \'id\' for an identity transformation for at leastone stack. The first stack with \'id\' transformation will be resampled as template.")
            ("transformations", value<string>(&parameters.frameRegfolder), "Use existing image-frame to volume transformations to initialize the reconstruction.")
            ("slice_transformations", value<string>(&parameters.sliceLocRegFolder), "Use existing slice-location transformations to initialize the reconstruction.")
            ("thickness", value<vector<double>>(&parameters.thickness)->multitoken(), "Give slice thickness. [Default: voxel size in z direction]")
            ("resolution", value<double>(&parameters.resolution), "Isotropic resolution of the volume. [Default: 0.75mm]")
            ("time_resolution", value<double>(&parameters.timeResolution), "Resolution on t-axis in milliseconds. [Default: average of stack resolution]")
            ("iterations", value<int>(&parameters.iterations), "Number of registration iterations [Default: 4]")
            ("sr_iterations", value<int>(&parameters.srIterations), "Number of SR reconstruction iterations for final iteration [Default: 21]")
            ("sigma", value<double>(&parameters.sigma), "Stdev for bias field. [Default: 12mm]")
            ("average", value<double>(&parameters.averageValue), "Average intensity value for stacks. [Default: 700]")
            ("delta", value<double>(&parameters.delta), "Parameter to define what is an edge. [Default: 150]")
            ("lambda", value<double>(&parameters.lambda), "Smoothing parameter. [Default: 0.02]")
            ("lastIter", value<double>(&parameters.lastIterLambda), "Smoothing parameter for last iteration. [Default: 0.01]")
            ("multires", value<int>(&parameters.levels), "Multiresolution smoothing with given number of levels. [Default: 3]")
            ("smooth_mask", value<double>(&parameters.smoothMask), "Smooth the mask to reduce artefacts of manual segmentation. [Default: 4mm]")
            ("force_exclude", value<vector<int>>(&parameters.forceExcludedSlices)->multitoken(), "Force exclusion of image-frames with these indices.")
            ("force_exclude_stack", value<vector<int>>(&parameters.forceExcludedStacks)->multitoken(), "Force exclusion of stacks with these indices.")
            ("no_stack_intensity_matching", "Switch off stack intensity matching.")
            ("stack_intensity_matching_only", "Switch off slice intensity matching.")
            ("no_intensity_matching", "Switch off all intensity matching.")
            ("no_robust_statistics", "Switch off robust statistics.")
            ("exclude_slices_only", bool_switch(&parameters.robustSlicesOnly), "Do not exclude individual voxels.")
            ("mask_slices_svr", "Mask slices for slice-to-volume registration only.")
            ("mask_slices_not_svr", "Mask slices everywhere except slice-to-volume registration.")
            ("mask_slices", "Mask slices always.")
            ("no_stack_zcrop", bool_switch(&parameters.noStackZCrop), "Disable cropping of stacks in the Z-direction to ROI i.e. keep all slices.")
            ("nmi_bins", value<int>(&parameters.nmiBins), "Number of NMI bins for registration. [Default: 16]")
            ("info", value<string>(&parameters.infoFilename), "File name for slice information in tab-separated columns.")
            ("log_prefix", value<string>(&parameters.logID), "Prefix for the log file.")
            ("save_slices", bool_switch(&parameters.saveSlicesFlag), "Save slices for future exclusion [Default: false]")
            ("save_evalulation_data", bool_switch(&parameters.saveEvaluationData), "Save outputs for use with evaluate-reconstruction tool. [Default: false]")
            ("debug", bool_switch(&parameters.debug), "Debug mode - save intermediate output results and writes debug information to log.")
            ("verbose", bool_switch(&parameters.verbose), "Verbose debugging mode - save highly-detailed itermediate output results and writes trace information to log.")
            ("profile", bool_switch(&parameters.profile), "Timings - output profiling timings");
            return opts;
    }

    void printUsage() {
        cout << "Usage: svr reconstruct [reconstructed] [stack_1] .. [stack_N] <options>\n" << endl;
        cout << "  [reconstructed]            Name for the reconstructed volume (Nifti format)" << endl;
        cout << "  [stack_1] .. [stack_N]     The input stacks (Nifti format)" << endl;
        ReconstructionParameters p;
        cout << getOptions(p) << endl;
    }

    int readInput(int argc, char **argv, ReconstructionParameters& parameters) {

        // Define required options
        options_description reqOpts;
        reqOpts.add_options()
            ("command", "reconstruct")
            ("reconstructed", value<string>(&parameters.outputName)->required(), "Name for the reconstructed volume (Nifti format)")
            ("stack", value<vector<string>>(&parameters.stackFiles)->multitoken()->required(), "The input stacks (Nifti format)");

        // Define positional options
        positional_options_description posOpts;
        posOpts.add("command", 1).add("reconstructed", 1).add("stack", -1);

        options_description opts = getOptions(parameters);
        options_description allOpts("Allowed options");
        allOpts.add(reqOpts).add(opts);

        // Parse arguments and catch errors
        variables_map vm;
        try {
            store(command_line_parser(argc, argv).options(allOpts).positional(posOpts)
                // Allow single dash (-) for long arguments
                .style(command_line_style::unix_style | command_line_style::allow_long_disguise).run(), vm);
            notify(vm);
        } catch (error& e) {
            // Delete -- from the argument name in the error message
            string err = e.what();
            size_t dashIndex = err.find("\'--");
            if (dashIndex != string::npos)
                err.erase(dashIndex + 1, 2);
            cerr << "Argument parsing error: " << err << "\n\n";
            printUsage();
            return 1;
        }
        
        if (parameters.verbose)
            logger::init_logging(logger::Level::TRACE, parameters.logID); 
        else if (parameters.debug) 
            logger::init_logging(logger::Level::DEBUG, parameters.logID);
        else 
            logger::init_logging(logger::Level::INFO, parameters.logID);

        // Read intensity matching flags and choose mode
        if      (vm.count("no_stack_intensity_matching")) { parameters.intensityMatchingMode = IM_MODE::FULL_RESCALE_EXCEPT_STACKS; }
        else if (vm.count("stack_intensity_matching_only")) { parameters.intensityMatchingMode = IM_MODE::RESCALE_STACKS_ONLY; }
        else if (vm.count("no_intensity_matching"))       { parameters.intensityMatchingMode = IM_MODE::NO_RESCALING; }
        else                                              { parameters.intensityMatchingMode = IM_MODE::FULL_RESCALING; }
        
        // Read slice masking flags and choose mode
        if      (vm.count("mask_slices_svr"))     { parameters.sliceMaskingMode = MASK_MODE::MASK_SLICES_ONLY_SVR; }
        else if (vm.count("mask_slices_not_svr")) { parameters.sliceMaskingMode = MASK_MODE::MASK_SLICES_NOT_SVR; } // 3D algorithm original behaivour.
        else if (vm.count("mask_slices"))         { parameters.sliceMaskingMode = MASK_MODE::MASK_SLICES_ALWAYS; }
        else                                      { parameters.sliceMaskingMode = MASK_MODE::NO_SLICE_MASKS; } //4D algorithm original behaviour

        // Switch off robust statistics
        if (vm.count("no_robust_statistics"))
            parameters.robustStatistics = false;

        if (!parameters.forceExcludedSlices.empty()) {
            stringstream ss;
            ss << parameters.forceExcludedSlices.size() << " force excluded slices: ";
            for (size_t i = 0; i < parameters.forceExcludedSlices.size(); i++)
                ss << parameters.forceExcludedSlices[i] << " ";
            logger::info(ss.str());
        }

        // Target stack
        if (vm.count("target_stack"))
            logger::info("Target stack no. is " + to_string(parameters.templateNumber) + " (zero-indexed stack no. " + to_string(--parameters.templateNumber) + ")");

        return 0;
    }

    int run(int argc, char **argv) {
        START_TIMING();

        ReconstructionParameters inputData;

        if (readInput(argc, argv, inputData) > 0)
            return 1;

        if (!inputData.Validate()) {
            logger::error("Invalid input parameters");
            return 1;
        }

        bool debug = inputData.debug;
        bool profile = inputData.profile;

        try {
            ReconstructionAlgorithm reconstructor(inputData);
            int result = reconstructor.Run();
            if (result > 0) {
                logger::error("Reconstruction failed!");
                return 1;
            }
        } catch(const std::invalid_argument& e) {
            string msg = "Invalid input parameters: ";
            msg.append(e.what());
            logger::error(msg);
            return 1;
        } catch (...) {
            logger::error("Reconstruction failed!");
            return 1;
        }

        END_TIMING("all");
        logger::info("End");
        logger::cleanup_logging();
        return 0;
    }

}