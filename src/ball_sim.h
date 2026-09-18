// src/ball_sim.h
#pragma once

#include "rolling_dynamics.h"
#include "table_kinematics.h"

#include <Eigen/Core>

namespace caliburn {

/// The ball the plate simulates: 40 mm across, 50 g, and rolling resistance of
/// one percent.
///
/// One definition, because the resistance is not decoration.  It is a Coulomb
/// term, so a plate tilted by less than `atan(0.01) = 0.573 deg` cannot start
/// the ball moving at all — the dead band a state feedback has no integral
/// term to escape.  A test that asserts where the loop stalls and a simulator
/// that stalls somewhere else would both be self-consistent and disagree.
inline constexpr BallParams kPlateBall{0.020, 0.050, 0.01};

/// Tilt angles as the rolling-ball model names them: `alpha` is the tilt that
/// drives the ball along y, `beta` the tilt that drives it along x.
struct BallTilt {
    double alpha;
    double beta;
};

/// Read the rolling model's tilt inputs off a table pose.
///
/// The two halves disagree, and this function is where they are reconciled:
/// `TableKinematics` uses R = Ry(theta) * Rx(phi), so the plate's local +Y edge
/// rises with phi while its local +X edge falls with theta.  The rolling model
/// accelerates the ball by +g*sin(alpha) along y and +g*sin(beta) along x.  So
/// theta carries straight through and phi has to flip.
BallTilt ballTiltFromPose(const TablePose& pose);

/// Advance the ball one RK4 step over a plate held at `pose`.
/// State is [x, y, vx, vy] in the plate's own frame.
///
/// `normal_accel` is the normal force per unit mass, N/m, in m/s^2, and it
/// scales the rolling resistance.  Held constant across the step, as the pose
/// is: both are what the plate is doing this frame.
Eigen::Vector4d stepBall(const RollingBallDynamics& dynamics,
                         const Eigen::Vector4d& state,
                         const TablePose& pose,
                         double normal_accel,
                         double dt);

/// The quasi-static overload: assumes the plate carries the ball's full weight.
/// Correct for a plate at rest, and the reading the tilt-convention tests want,
/// since they ask what a held pose does and not what a moving plate does.
Eigen::Vector4d stepBall(const RollingBallDynamics& dynamics,
                         const Eigen::Vector4d& state,
                         const TablePose& pose,
                         double dt);

/// True while the ball's contact patch is still over a disc plate.
/// The plate is a disc, not the square that `RollingBallDynamics::on_plate`
/// assumes, so the bound is radial and inset by the ball's own radius.
bool ballOnPlate(const Eigen::Vector4d& state,
                 double plate_radius,
                 double ball_radius);

}  // namespace caliburn
