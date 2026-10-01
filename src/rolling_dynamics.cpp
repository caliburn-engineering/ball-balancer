// From caliburn reference/simulation/rolling_dynamics.cpp @ ee2cf6e
#include "rolling_dynamics.h"

#include <cmath>

namespace caliburn {

RollingBallDynamics::RollingBallDynamics(const BallParams& ball, const PlateParams& plate)
    : ball_(ball), plate_(plate) {}

Eigen::Vector4d RollingBallDynamics::derivatives(const Eigen::Vector4d& state,
                                                  double alpha, double beta) const {
    return derivatives(state, alpha, beta, plate_.gravity);
}

Eigen::Vector4d RollingBallDynamics::derivatives(const Eigen::Vector4d& state,
                                                  double alpha, double beta,
                                                  double normal_accel) const {
    double vx = state(2);
    double vy = state(3);
    double g = plate_.gravity;
    double rf = rolling_factor();
    double crr = ball_.rolling_friction;

    // The surface can push, never pull: a negative normal force is not a
    // backwards friction, it is a ball that is no longer touching. Callers that
    // model separation decide that before they get here; clamping keeps a
    // caller that does not from being handed a friction that accelerates.
    double n = (normal_accel > 0.0) ? normal_accel : 0.0;

    // Gravity component along plate surface
    double grav_x = g * std::sin(beta);
    double grav_y = g * std::sin(alpha);

    // Rolling friction (opposes velocity), scaled by the normal force
    double fric_x = (vx != 0.0) ? crr * n * std::copysign(1.0, vx) : 0.0;
    double fric_y = (vy != 0.0) ? crr * n * std::copysign(1.0, vy) : 0.0;

    // Effective acceleration with rolling factor
    double ax = rf * (grav_x - fric_x);
    double ay = rf * (grav_y - fric_y);

    return Eigen::Vector4d(vx, vy, ax, ay);
}

bool RollingBallDynamics::on_plate(const Eigen::Vector4d& state) const {
    return std::abs(state(0)) < plate_.half_width &&
           std::abs(state(1)) < plate_.half_width;
}

const BallParams& RollingBallDynamics::ball_params() const { return ball_; }
const PlateParams& RollingBallDynamics::plate_params() const { return plate_; }

}  // namespace caliburn
