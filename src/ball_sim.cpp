// src/ball_sim.cpp
#include "ball_sim.h"

#include "rk4.h"

#include <cmath>

namespace caliburn {

BallTilt ballTiltFromPose(const TablePose& pose) {
    return BallTilt{-pose.phi, pose.theta};
}

Eigen::Vector4d stepBall(const RollingBallDynamics& dynamics,
                         const Eigen::Vector4d& state,
                         const TablePose& pose,
                         double normal_accel,
                         double dt) {
    const BallTilt tilt = ballTiltFromPose(pose);
    DerivativeFn f = [&](double, const Eigen::VectorXd& y) -> Eigen::VectorXd {
        return dynamics.derivatives(Eigen::Vector4d(y), tilt.alpha, tilt.beta,
                                    normal_accel);
    };
    Eigen::VectorXd next = rk4_step(Eigen::VectorXd(state), 0.0, dt, f);
    return Eigen::Vector4d(next);
}

Eigen::Vector4d stepBall(const RollingBallDynamics& dynamics,
                         const Eigen::Vector4d& state,
                         const TablePose& pose,
                         double dt) {
    return stepBall(dynamics, state, pose,
                    dynamics.plate_params().gravity, dt);
}

bool ballOnPlate(const Eigen::Vector4d& state,
                 double plate_radius,
                 double ball_radius) {
    const double r = std::hypot(state(0), state(1));
    return r <= plate_radius - ball_radius;
}

}  // namespace caliburn
