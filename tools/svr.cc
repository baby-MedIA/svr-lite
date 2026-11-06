// Copyright © King's College London
// SPDX-License-Identifier: Apache-2.0
//
// This file has been modified from the original source.

// SVR

#define SVR_TOOL
#include "svr/Profiling.h"
#include "reconstruct.h"
#include "evaluate.h"

using namespace std;
using namespace svr;
using namespace boost::program_options;

options_description getOptions() {
    options_description opts("Options");
    opts.add_options()
        ("help", "Usage help with svr commands.");
    return opts;
}

void printUsage() {
    cout << "Usage: svr <command>" << endl;
    cout << "Commands:" << endl;
    cout << "   reconstruct           Create a slice-to-volume reconstruction." << endl;
    cout << "   evaluate              Evaluate the reconstructed volume against specified stacks." << endl;
    cout << getOptions() << endl;
    cout << "Use --help <command> for further usage help." << endl;
}

int main(int argc, char **argv) {
    options_description opts = getOptions();
    variables_map vm;

    options_description reqOpts;
    string command;
    reqOpts.add_options()
        ("command", value<string>(&command), "SVR command to execute.");
        
    options_description allOpts("Allowed options");
    allOpts.add(reqOpts).add(opts);
    
    positional_options_description posOpts;
    posOpts.add("command", 1);

    // only view first command at this level
    int svrArgC = (argc > 2) ? 2 : argc;

    try {
        store(command_line_parser(svrArgC, argv).options(allOpts).positional(posOpts)
            .style(command_line_style::unix_style | command_line_style::allow_long_disguise).run(),vm);
        notify(vm);
    } catch (error& e) {
        // Delete -- from the argument name in the error message
        string err = e.what();
        cerr << "Argument parsing error: " << err << "\n\n";
        printUsage();
        return 1;
    }

    // convert command to upper case for usability
    transform(command.begin(), command.end(), command.begin(), ::toupper);

    // Run help
    if (vm.count("help")) {
        int i = 0;
        int count = 0;

        // read all args and print appropriate usage message
        while (i < argc) {
            string arg = argv[i];
            transform(arg.begin(), arg.end(), arg.begin(), ::toupper);
            if (arg == "RECONSTRUCT")   { reconstruct_tool::printUsage(); ++count; }
            if (arg == "EVALUATE")      { ++count; } //evaluate_tool::printUsage(); ++count; }
            ++i;
        }

        // give general usage help if no other help provided
        if (count <= 0){
            printUsage();
        }
    }

    // Run tools
    else if (command=="RECONSTRUCT")    {   reconstruct_tool::run(argc, argv);  }
    else if (command=="EVALUATE")       {   evaluate_tool::run(argc, argv);     }

    // default: usage
    else {
        cout << "print usage" << endl;
        printUsage();
    }
}