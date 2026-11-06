// Copyright © King's College London
// SPDX-License-Identifier: Apache-2.0
//
// This file has been modified from the original source.

#pragma once

// SVR
#include "svr/ReconstructionModel.h"

using namespace mirtk;
using namespace svr;
using namespace svr::Utility;

namespace svr::Parallel {
        // class for SStep (local structure-based outlier rejection)

    class SStep {
        ReconstructionModel *reconstructor;
        size_t N;

    public:
        SStep(ReconstructionModel *reconstructor, size_t N) : reconstructor(reconstructor),
            N(N) {}

        void operator()(const blocked_range<size_t>& r) const {
            for (size_t inputIndex = r.begin(); inputIndex != r.end(); inputIndex++) {

                reconstructor->SliceSSIMMap(inputIndex);

            } //end of loop for a slice inputIndex
        }

        void operator()() {
            parallel_for(blocked_range<size_t>(0, N), *this);
        }
    };

}