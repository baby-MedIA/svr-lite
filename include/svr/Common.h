// Copyright © King's College London
// SPDX-License-Identifier: Apache-2.0
//
// This file has been modified from the original source.

#pragma once

// MIRTK
#include "mirtk/Common.h"
#include "mirtk/Options.h"
#include "mirtk/IOConfig.h"
#include "mirtk/Array.h"
#include "mirtk/Arith.h"
#include "mirtk/Math.h"
#include "mirtk/Point.h"
#include "mirtk/BaseImage.h"
#include "mirtk/GenericImage.h"
#include "mirtk/ImageReader.h"
#include "mirtk/Dilation.h"
#include "mirtk/Erosion.h"
#include "mirtk/GaussianBlurring.h"
#include "mirtk/Resampling.h"
#include "mirtk/ResamplingWithPadding.h"
#include "mirtk/LinearInterpolateImageFunction.hxx"
#include "mirtk/GaussianBlurringWithPadding.h"
#include "mirtk/GenericRegistrationFilter.h"
#include "mirtk/Transformation.h"
#include "mirtk/HomogeneousTransformation.h"
#include "mirtk/RigidTransformation.h"
#include "mirtk/ImageTransformation.h"
#include "mirtk/MultiLevelFreeFormTransformation.h"
#include "mirtk/FreeFormTransformation.h"
#include "mirtk/LinearFreeFormTransformation3D.h"
#include "mirtk/BSplineFreeFormTransformationSV.h"

// C++ Standard
#include <algorithm>
#include <pthread.h>
#include <queue>
#include <set>
#include <thread>

// Boost
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/program_options.hpp>
#include <boost/timer/timer.hpp>

// OpenMP
#include <omp.h>

// SVR
#include "svr/MeanShift.h"
#include "svr/Utility.h"
#include "svr/Logger.h"

namespace svr {
    struct POINT3D {
        short x;
        short y;
        short z;
        double value;
    };

    enum RECON_TYPE { _3D, _1D, _interpolate };

    typedef Array<POINT3D> VOXELCOEFFS;
    typedef Array<Array<VOXELCOEFFS>> SLICECOEFFS;

    /// PI
    constexpr double PI = 3.14159265358979323846;
}
