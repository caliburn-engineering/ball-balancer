// src/analysis/loop_locus.h
#pragma once

#include "linear_system.h"
#include "controller_builders.h"
#include "pole_zero.h"

#include <vector>

namespace caliburn {

// Root locus of 1 + kappa*L(s) = 0, where kappa scales the WHOLE compensator on
// rule C's channel (input_j, output_i) — every loop with out == output_i and
// in == input_j, so duplicate pairings sum naturally and no loop selector is
// needed.  Other loops are held at their current tuning.
//
// This is a true root locus, not a named-gain parameter sweep: with Ki != 0 a
// Kp sweep does not open the loop at zero, branches do not start at open-loop
// poles, and the asymptote rules do not hold.  Scaling the whole compensator
// keeps all of that meaning what the help text says.
//
// This module sits ABOVE pole_zero.h so that module stays what it is today: a
// general routine over any LinearSystem with no knowledge of controllers.
//
// kappa is log-spaced over [kappa_min, kappa_max], kappa_min > 0.  Returns
// empty when the channel carries no loop.
std::vector<RootLocusPoint> computeLoopLocus(
    const LinearSystem& plant,
    const std::vector<Loop>& loops,
    LoopKind kind,
    int output_i, int input_j,
    double kappa_min, double kappa_max, int num_points);

// The selected loop's gain margin with every other loop held at its tuning:
// the first sign change in max Re(pole), linearly interpolated.  Returns a
// negative value when there is no crossing in range — never a fabricated number.
double loopGainMargin(const std::vector<RootLocusPoint>& locus);

}  // namespace caliburn
