// tests/test_ball_sim.cpp
//
// The ball's tilt convention is the one thing in the merged plate view that
// cannot be checked by reading the code: RollingBallDynamics takes (alpha,
// beta), TableKinematics produces (phi, theta), and the two are neither in the
// same order nor the same sign.  A wrong mapping still compiles, still moves
// the ball, and only shows up as a ball that rolls uphill.
#include "ball_sim.h"
#include "test_helpers.h"

#include <cmath>

using namespace caliburn;

namespace {

RollingBallDynamics makeDynamics(double rolling_friction = 0.0) {
    BallParams ball{0.02, 0.05, rolling_friction};
    PlateParams plate{0.30, 9.81};
    return RollingBallDynamics(ball, plate);
}

constexpr double kRollingFactor = 5.0 / 7.0;

// --- Pitch (theta, about +Y) drops the +X edge, so the ball runs to +x ---
void test_pitch_rolls_ball_towards_positive_x() {
    const auto dyn = makeDynamics();
    const double theta = 0.10;
    const double dt = 0.02;

    TablePose pose{0.0, theta, 0.12};
    Eigen::Vector4d x = stepBall(dyn, Eigen::Vector4d::Zero(), pose, dt);

    // RK4 is exact for a constant acceleration, so this is an equality, not a
    // trend: a = (5/7) g sin(theta).
    ASSERT_NEAR(x(2), kRollingFactor * 9.81 * std::sin(theta) * dt, 1e-12);
    ASSERT_NEAR(x(3), 0.0, 1e-12);
    ASSERT_TRUE(x(0) > 0.0);
    ASSERT_NEAR(x(1), 0.0, 1e-12);
}

// --- Roll (phi, about +X) lifts the +Y edge, so the ball runs to -y ---
void test_roll_rolls_ball_towards_negative_y() {
    const auto dyn = makeDynamics();
    const double phi = 0.10;
    const double dt = 0.02;

    TablePose pose{phi, 0.0, 0.12};
    Eigen::Vector4d x = stepBall(dyn, Eigen::Vector4d::Zero(), pose, dt);

    ASSERT_NEAR(x(3), -kRollingFactor * 9.81 * std::sin(phi) * dt, 1e-12);
    ASSERT_NEAR(x(2), 0.0, 1e-12);
    ASSERT_TRUE(x(1) < 0.0);
    ASSERT_NEAR(x(0), 0.0, 1e-12);
}

// --- The mapping itself, stated once ---
void test_tilt_mapping() {
    BallTilt t = ballTiltFromPose(TablePose{0.3, -0.2, 0.12});
    ASSERT_NEAR(t.alpha, -0.3, 1e-15);
    ASSERT_NEAR(t.beta, -0.2, 1e-15);
}

// --- A level plate leaves a resting ball alone ---
void test_level_plate_is_an_equilibrium() {
    const auto dyn = makeDynamics();
    Eigen::Vector4d x(0.05, -0.03, 0.0, 0.0);
    for (int i = 0; i < 100; ++i)
        x = stepBall(dyn, x, TablePose{0.0, 0.0, 0.12}, 0.01);

    ASSERT_NEAR(x(0), 0.05, 1e-12);
    ASSERT_NEAR(x(1), -0.03, 1e-12);
    ASSERT_NEAR(x(2), 0.0, 1e-12);
    ASSERT_NEAR(x(3), 0.0, 1e-12);
}

// --- The plate is a disc, so the bound is radial, not per-axis ---
void test_on_plate_is_radial() {
    const double R = 0.30, r = 0.02;
    ASSERT_TRUE(ballOnPlate(Eigen::Vector4d(0.0, 0.0, 0, 0), R, r));
    ASSERT_TRUE(ballOnPlate(Eigen::Vector4d(0.27, 0.0, 0, 0), R, r));

    // Inside the square |x|,|y| < R, but outside the disc.
    ASSERT_TRUE(!ballOnPlate(Eigen::Vector4d(0.25, 0.25, 0, 0), R, r));

    // The ball rests on its contact patch, so it leaves at R - r, not R.
    ASSERT_TRUE(!ballOnPlate(Eigen::Vector4d(0.29, 0.0, 0, 0), R, r));
}

}  // namespace

int main() {
    test_pitch_rolls_ball_towards_positive_x();
    test_roll_rolls_ball_towards_negative_y();
    test_tilt_mapping();
    test_level_plate_is_an_equilibrium();
    test_on_plate_is_radial();
    std::printf("test_ball_sim: all passed\n");
    return 0;
}
