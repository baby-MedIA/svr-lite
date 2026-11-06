// Copyright © King's College London
// SPDX-License-Identifier: Apache-2.0
//
// This file has been modified from the original source.

// SVR
#include "svr/ReconstructionAlgorithm.h"
#include "svr/Evaluation.h"

 
#define SVR_TOOL
#include "svr/Profiling.h"

using namespace std;
using namespace mirtk;
using namespace svr;
using namespace boost::program_options;

namespace evaluate_tool {

    void printUsage(options_description opts) {
        cout << "Usage: svr evaluate [reconstructed image] [Stack to evaluate] <options>\n" << endl;
        cout << "  [reconstructed]            Reconstructed image (Nifti format)" << endl;
        cout << "  [stack] .. [stack_N]       Stacks to evaluate (Nifti format)" << endl;
        cout << opts << endl;
    }

    int run(int argc, char **argv) {
        // START_TIMING();

        // Define input variables
        vector<string> stackRegFolder;
        double thickness = -1;
        string stackFile, maskFile, reconVolume, sliceLocRegFolder;
        bool stackIntensityMatching = false;
        bool debug = false;
        bool verbose = false;
        bool profile = false;
        int iterations = -1;

        // Define optional parameters
        options_description opts("Options");
        opts.add_options()
            ("thickness", value<double>(&thickness), "Give slice thickness. [Default: voxel size in z direction]")
            ("iterations", value<int>(&iterations), "Number of registration iterations [Default: Auto-detect convergence]")
            ("stack_intensity_matching", bool_switch(&stackIntensityMatching), "Perform stack intensity matching. [Default: off]")
            ("debug", bool_switch(&debug), "Debug mode - save intermediate output results and writes debug information to log. [Default: false]")
            ("verbose", bool_switch(&verbose), "Verbose debugging mode - save highly-detailed itermediate output results and writes trace information to log. [Default: false]")
            ("profile", bool_switch(&profile), "Timings - output profiling timings. [Default: false]");

        // Define required options
        options_description reqOpts;
        reqOpts.add_options()
            ("command", "evaluate")
            ("reconstructed", value<string>(&reconVolume)->required(), "Name for the reconstructed volume (Nifti format)")
            ("stack", value<string>(&stackFile)->required(), "Input stack (Nifti format)")
            ("mask", value<string>(&maskFile)->required(), "Mask that defines evaluation region of interest (may be smaller than reconstructed ROI).");

        // Define positional options
        positional_options_description posOpts;
        posOpts.add("command", 1).add("reconstructed", 1).add("stack", 1).add("mask", 1);

        //options_description opts = getOptions(parameters);
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
            printUsage(opts);
            return 1;
        }

        // Run evaluation
        Evaluation eval(debug, verbose, profile);
        eval.run(stackFile, thickness, maskFile, reconVolume, stackIntensityMatching, iterations);

        return 0;
    }

}