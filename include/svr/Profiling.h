// Copyright © King's College London
// SPDX-License-Identifier: Apache-2.0
//
// This file has been modified from the original source.

#pragma once

#include <iostream>
#include <iomanip>
#include <ctime>
#include <chrono>
#include <boost/timer/timer.hpp>
#include "svr/Logger.h"

#define START_TIMING() boost::timer::cpu_timer timer; timer.start();

#define RESET_TIMING() timer.start();

#define OUTPUT_TIMING(section)\
    do {\
        timer.stop();\
        std::stringstream t_out;\
        t_out << "Time for " << section << ": " << timer.format(3,"%w secs ");\
        logger::info(t_out.str());\
        std::stringstream t_prof;\
        t_prof << section << "," << timer.format(3,"%w") << ",secs";\
        logger::profile(t_prof.str());\
    }while(false);


#ifdef  SVR_TOOL
    #define END_TIMING(section) if (profile) { OUTPUT_TIMING(section) }
#else
    #define END_TIMING(section) if (_profile) { OUTPUT_TIMING(section) }
#endif
