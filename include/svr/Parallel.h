// Copyright © King's College London
// SPDX-License-Identifier: Apache-2.0
//
// This file has been modified from the original source.

#pragma once

/// Include parallel functors
#include "svr/parallel/Average.h"
#include "svr/parallel/AdaptiveRegularisation.h"
#include "svr/parallel/Bias.h"
#include "svr/parallel/CalculateCorrectedSlices.h"
#include "svr/parallel/CalculateError.h"
#include "svr/parallel/CoeffInit.h"
#include "svr/parallel/EStep.h"
#include "svr/parallel/MStep.h"
#include "svr/parallel/SStep.h"
#include "svr/parallel/NormaliseBias.h"
#include "svr/parallel/Scale.h"
#include "svr/parallel/SimulateSlices.h"
#include "svr/parallel/SimulateStacks.h"
#include "svr/parallel/SliceToVolumeRegistration.h"
#include "svr/parallel/StackRegistrations.h"
#include "svr/parallel/SuperResolution.h"
