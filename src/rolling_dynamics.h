// From caliburn reference/simulation/rolling_dynamics.h @ ee2cf6e
#pragma once

#include <Eigen/Core>

namespace caliburn {

struct BallParams {
    double radius;           // ball radius [m]
    double mass;             // ball mass [kg]
    double rolling_friction; // coefficient of rolling resistance
};

struct PlateParams {
    double half_width;   // plate half-width [m] (square plate)
    double gravity;      // gravitational acceleration [m/s²]
};

class RollingBallDynamics {
public:
    RollingBallDynamics(const BallParams& ball, const PlateParams& plate);

    // Compute state derivative: state = [x, y, vx, vy], tilt = [alpha, beta] in radians.
    // Returns [vx, vy, ax, ay].
    //
    // `normal_accel` is the normal force per unit mass, N/m, in m/s^2. Rolling
    // resistance is a normal-force effect — the contact patch deforms in
    // proportion to how hard the surface is pressed, not to how strong gravity
    // is — so the retarding acceleration is c_rr * (N/m), not c_rr * g. The two
    // agree only on a level plate at rest, which is the case this model was
    // first written for and the assumption that was buried in it.
    //
    // A plate that is tilted, heaving or rotating presses the ball with
    // something other than its weight, and a plate that is falling away from
    // the ball presses it with almost nothing. See issue #23.
    Eigen::Vector4d derivatives(const Eigen::Vector4d& state,
                                double alpha, double beta,
                                double normal_accel) const;

    // The quasi-static reading: assumes the surface carries the ball's full
    // weight, N/m = g. Correct for a level plate at rest and the right choice
    // for linearising about that equilibrium; wrong for a plate in motion.
    Eigen::Vector4d derivatives(const Eigen::Vector4d& state,
                                double alpha, double beta) const;

    // Check if ball is on plate.
    bool on_plate(const Eigen::Vector4d& state) const;

    // The 5/7 effective acceleration factor for a uniform solid sphere.
    static constexpr double rolling_factor() { return 5.0 / 7.0; }

    const BallParams& ball_params() const;
    const PlateParams& plate_params() const;

private:
    BallParams ball_;
    PlateParams plate_;
};

}  // namespace caliburn
