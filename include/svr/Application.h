// Copyright © King's College London
// SPDX-License-Identifier: Apache-2.0
//
// This file has been modified from the original source.

#pragma once

// SVR
#include <filesystem>
#include "svr/Logger.h"

namespace svr {

    /***
     * Static base class for SVR application
     * Defines common method interface and variables required by logger and profiler
     */
    class Application {
        protected:
            // Debug mode (user debugging)
            bool _debug = false;

            // Verbose debugging mode (developer debugging)
            bool _verbose = false;

            // Profilling mode
            bool _profile = false;

            // Set up logger to work with debug/verbose mode flags
            inline void ConfigureLogging(string logID) {
                if (_verbose)
                    logger::init_logging(logger::Level::TRACE, logID);
                else if (_debug)
                    logger::init_logging(logger::Level::DEBUG, logID);
                else
                    logger::init_logging(logger::Level::INFO, logID);
            }
    };
}