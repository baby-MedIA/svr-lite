// Copyright © King's College London
// SPDX-License-Identifier: Apache-2.0
//
// This file has been modified from the original source.

#pragma once

#include <iostream>
#include "spdlog/spdlog.h"
#include "spdlog/sinks/basic_file_sink.h"

using namespace std;

namespace logger     {
    static const string LOGNAME = "svr";
    static const string LOGFILE = "main.log";

    static const string PROFILELOG = "profiling";
    static const string PROFILELOGFILE = "profile.log";

    // Add further log levels as need requires
    enum Level { TRACE, DEBUG, INFO, ERROR }; 

    static const Level OUTPUT_AT = Level::INFO; // Output messages to cout if this level or higher.

    inline void init_logging(Level logLevel, string prefix="") {
        if (spdlog::get(LOGNAME) == nullptr) {
            try 
            {
                auto logger = spdlog::basic_logger_mt(LOGNAME, prefix+LOGFILE);
                logger->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%t] [%l] %v");

                switch (logLevel) {
                    case ERROR: logger->set_level(spdlog::level::err);   break;
                    case INFO:  logger->set_level(spdlog::level::info);  break;
                    case DEBUG: logger->set_level(spdlog::level::debug); break;
                    case TRACE: logger->set_level(spdlog::level::trace); break;
                    default:    logger->set_level(spdlog::level::info); break; 
                }
            }
            catch (const spdlog::spdlog_ex &ex) {
                std::cerr << "Log initialisation failed: " << ex.what() << endl;
                throw ex;
            }
        }
    }

    inline void trace(string msg) {
        auto logger = spdlog::get(LOGNAME);
        if (logger != nullptr) 
            logger->trace(msg);
        if (Level::TRACE >= OUTPUT_AT)
            cout << msg << endl;
    }

    inline void debug(string msg) {
        auto logger = spdlog::get(LOGNAME);
        if (logger != nullptr) 
            logger->debug(msg);

        if (Level::DEBUG >= OUTPUT_AT)
            cout << msg << endl; 
    }

    inline void info(string msg) {
        auto logger = spdlog::get(LOGNAME);
        if (logger != nullptr) 
            logger->info(msg);

        if (Level::INFO >= OUTPUT_AT)
            cout << msg << endl;
    }

    inline void error(string msg) {
        auto logger = spdlog::get(LOGNAME);
        if (logger != nullptr) 
            logger->error(msg);

        if (Level::ERROR >= OUTPUT_AT)
            cerr << msg << endl; 
    } 

    inline void cleanup_logging(){
        spdlog::drop_all();
    }

    inline void flush() {
        for (auto name : {LOGNAME, PROFILELOG}) {
            auto logger = spdlog::get(name);
            if (logger != nullptr) 
                logger->flush();
        }
    }

    /***
     * Special log functions for profiling:
     * Set up profile log.
     */
    inline void init_profile_log() {
        if (spdlog::get(PROFILELOG) == nullptr) {
            try 
            {
                auto logger = spdlog::basic_logger_mt(PROFILELOG, PROFILELOGFILE);
                logger->set_pattern("%v");
                logger->set_level(spdlog::level::info);  
            }
            catch (const spdlog::spdlog_ex &ex) {
                std::cerr << "Profile log initialisation failed: " << ex.what() << endl;
                throw ex;
            }
        }
    }

    /***
     * Special log functions for profiling:
     * Write to profile log.
     */
    inline void profile(string msg) {
        auto logger = spdlog::get(PROFILELOG);
        if (logger != nullptr) 
            logger->info(msg);
    }
}