#pragma once

#include "linear_system.h"
#include "linearizer.h"

namespace caliburn {

struct BallParams;
struct PlateParams;

/// Hand-derived linear model for ball rolling on tilted plate.
/// Operating point: ball at centre, plate level.
/// State: [x, y, vx, vy], Input: [alpha, beta]
/// Output: [x, y] (position only)
LinearSystem ballBalancerLinearModel(double gravity = 9.81);

/// Wrap RollingBallDynamics::derivatives as a NonlinearFn for comparison.
/// Drops the time argument (linearization is time-invariant).
NonlinearFn ballBalancerNonlinearFn(double ball_radius = 0.02,
                                    double ball_mass = 0.05,
                                    double rolling_friction = 0.0,
                                    double plate_half_width = 0.15,
                                    double gravity = 9.81);

}  // namespace caliburn
