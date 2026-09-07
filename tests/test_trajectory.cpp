// tests/test_trajectory.cpp
//
// The ball following a path, against the nonlinear plate.  `test_setpoint_path`
// pins the shapes; this pins what the controller does with them, which is the
// part a visitor actually watches.
//
// Every combination the UI can be put into is swept, because "the ball stays
// on the plate at every offered speed" is a claim about the SLIDERS, and a
// slider setting that loses the ball is not a choice on offer, it is a trap.
// See issue #24.
#include "setpoint_path.h"

#include "auto_balance.h"
#include "ball_contact.h"
#include "ball_sim.h"
#include "cascade_fixture.h"
#include "table_kinematics.h"
#include "test_helpers.h"

#include <cmath>

using namespace caliburn;

namespace {

constexpr double kBallRadius = kFixtureBallRadius;

struct Track {
    double max_radius = 0.0;   ///< [m] furthest the BALL got from centre
    double mean_error = 0.0;   ///< [m] after the first lap
    double max_error = 0.0;
    bool lost = false;
};

/// PlateView::step's calls in PlateView::step's order, with the setpoint driven
/// by a path and the loop given the path's velocity as its reference velocity.
///
/// `feedforward` chooses whether that velocity is SUPPLIED, and no longer how
/// it is applied: the arithmetic lives in `legCommand`, where the application
/// gets it from too.  This harness used to carry its own copy — a second
/// implementation of the loop is what let a real bug hide once already (#23),
/// and a harness that computes the thing it is meant to be measuring is not
/// evidence about the application.
Track follow(const ModelEntry& e, const Eigen::MatrixXd& K,
             const SetpointPath& path, bool feedforward, double duration) {
    const TableParams tp = cascadeMechanism(e.params);
    const TableKinematics tk(tp);
    const double g = cascadeGravity(e.params);
    const RollingBallDynamics dynamics(kPlateBall, PlateParams{tp.R_table, g});

    AutoBalanceDesign d;
    d.K = K;
    d.home_leg_rad = cascadeHomeLegAngle(e.params);
    d.servo_tau = cascadeServoTau(e.params);
    d.alpha_min_rad = tp.alpha_min;
    d.alpha_max_rad = tp.alpha_max;

    std::array<double, 3> alpha = {d.home_leg_rad, d.home_leg_rad, d.home_leg_rad};
    TablePose pose = tk.home_pose(d.home_leg_rad);
    PlateMotion pm = plateMotion(tk, pose, alpha, {0.0, 0.0, 0.0}), pm_prev = pm;

    BallState ball;
    const Eigen::Vector2d start = pathPoint(path, 0.0);
    ball.rolling << start(0), start(1), 0.0, 0.0;   // starts on the path

    const double dt = 1.0 / 60.0;
    Track r;
    double t = 0.0, phase = 0.0, sum = 0.0;
    int n = 0;
    for (int k = 0; k < static_cast<int>(duration / dt); ++k) {
        // `stepPath`, which is what the application calls and in the order it
        // calls it: read at the phase the frame opened on, advance afterwards.
        // Not re-derived here — a harness that recomputes the thing it is
        // measuring stops being evidence about the application.
        const PathStep step = stepPath(path, phase, dt);
        const Eigen::Vector2d sp = step.point;
        const Eigen::Vector2d spv = step.velocity;
        phase = step.next_phase;

        const Eigen::Matrix<double, 6, 1> bp = plateFrame(ball, pm, kBallRadius);
        Eigen::Vector4d seen(bp(0), bp(1), bp(3), bp(4));
        if (ball.airborne) {
            const Eigen::Vector2d land = predictedLanding(bp, kBallRadius, g);
            seen << land(0), land(1), 0.0, 0.0;
        }
        BallReference ref;
        ref.position = sp;
        if (feedforward && !ball.airborne) ref.velocity = spv;
        const LegCommand c = legCommand(tk, d, alpha, seen, ref);
        std::array<double, 3> adot{};
        for (int j = 0; j < 3; ++j)
            adot[j] = (c.alpha_rad[j] - alpha[j]) / d.servo_tau;
        alpha = stepServosOnPlate(tk, alpha, c.alpha_rad, d.servo_tau, dt);
        const FKResult fk = tk.solve_pose(alpha, pose);
        if (fk.converged) pose = fk.pose;
        pm_prev = pm;
        pm = plateMotion(tk, pose, alpha, adot,
                         servoAccel(adot, d.servo_tau), &pm_prev, dt);
        ball = stepBallContact(dynamics, ball, pm, pm_prev, pose, kBallRadius, g, dt);
        t += dt;

        const Eigen::Matrix<double, 6, 1> now = plateFrame(ball, pm, kBallRadius);
        r.max_radius = std::max(r.max_radius, std::hypot(now(0), now(1)));
        if (t > path.period_s) {          // after one lap, so the start is not counted
            const double err = std::hypot(now(0) - sp(0), now(1) - sp(1));
            r.max_error = std::max(r.max_error, err);
            sum += err;
            ++n;
        }
        if (!ballOnPlate(Eigen::Vector4d(now(0), now(1), now(3), now(4)),
                         tp.R_table, kBallRadius)) {
            r.lost = true;
            return r;
        }
    }
    r.mean_error = n ? sum / n : 0.0;
    return r;
}

// The acceptance criterion, over the whole envelope the sliders offer: every
// shape, the size range end to end, and at each size the fastest lap the UI
// will allow as well as the slowest.
void test_every_offered_setting_keeps_the_ball() {
    const auto models = getBuiltinModels();
    const auto& e = cascadeModel(models);
    const Eigen::MatrixXd K = defaultGain(e);
    const TableParams tp = cascadeMechanism(e.params);
    const double rim = tp.R_table - kBallRadius;

    double worst_radius = 0.0, worst_error = 0.0;
    for (PathShape s : {PathShape::Circle, PathShape::Square, PathShape::Triangle}) {
        for (double r_mm : {20.0, 60.0, 120.0, kMaxPathRadius * 1000.0}) {
            SetpointPath p;
            p.shape = s;
            p.radius_m = r_mm * 1e-3;
            // The fastest the UI allows at this size, and the slowest.
            for (double T : {minPeriod(p), 30.0}) {
                p.period_s = T;
                const Track t = follow(e, K, p, true, 2.0 * T + 4.0);
                ASSERT_TRUE(!t.lost);
                // And not merely on the plate — comfortably on it.  The ball
                // swings wide of a corner, and that overshoot is what the
                // radius cap exists to leave room for.
                ASSERT_TRUE(t.max_radius < 0.8 * rim);
                worst_radius = std::max(worst_radius, t.max_radius);
                worst_error = std::max(worst_error, t.mean_error);
            }
        }
    }
    // Measured: the ball reaches 195 mm on the largest path, against the 280 mm
    // the plate has, and its mean error is 36 mm at the very fastest setting
    // the sliders offer — which is the setting that is SUPPOSED to look hard.
    ASSERT_TRUE(worst_radius < 0.24);
    ASSERT_TRUE(worst_error < 0.05);
}

// Velocity feedforward, as the measurement that decided it.  The reference
// state claims the ball should be at the setpoint AND stationary, which is
// false the moment the setpoint moves; without the feedforward the loop spends
// its effort fighting the motion it was asked for, and the ball hangs back
// near the centre while the setpoint goes round without it.
void test_feedforward_is_worth_having() {
    const auto models = getBuiltinModels();
    const auto& e = cascadeModel(models);
    const Eigen::MatrixXd K = defaultGain(e);

    SetpointPath p;
    p.shape = PathShape::Circle;
    p.radius_m = 0.12;
    p.period_s = 10.0;

    const Track with = follow(e, K, p, true, 24.0);
    const Track without = follow(e, K, p, false, 24.0);

    ASSERT_TRUE(!with.lost && !without.lost);
    // Measured 3.66 mm against 35.52 mm — nearly ten times.  The assertion is
    // five, so that a change which merely halves the benefit still fails.
    ASSERT_TRUE(with.mean_error * 5.0 < without.mean_error);
    // And without it the error is a third of the radius: the ball is not
    // following the circle so much as sitting inside it.
    ASSERT_TRUE(without.mean_error > 0.25 * p.radius_m);
}

// The corner is the point of the cornered shapes.  A square is tracked worse
// than a circle of the same size and lap time, and it has to be — the corner
// is a step in the reference velocity, and a bandwidth-limited loop rounds it.
void test_a_corner_is_harder_than_a_curve() {
    const auto models = getBuiltinModels();
    const auto& e = cascadeModel(models);
    const Eigen::MatrixXd K = defaultGain(e);

    SetpointPath circle;
    circle.shape = PathShape::Circle;
    circle.radius_m = 0.12;
    circle.period_s = 8.0;
    SetpointPath square = circle;
    square.shape = PathShape::Square;

    const Track c = follow(e, K, circle, true, 20.0);
    const Track q = follow(e, K, square, true, 20.0);
    ASSERT_TRUE(!c.lost && !q.lost);
    // The square's worst error is at its corners: measured 14.2 mm against the
    // circle's 5.1 mm, on the same size and the same lap time.  Nearly three
    // times, and it is the corner that does it.
    //
    // The assertion is a RATIO rather than either number, because the square's
    // figure is a peak sampled at 60 Hz right where the reference velocity
    // steps: which frame falls nearest the corner decides it, and a phase
    // perturbation at the last bit of a double swings it between 14.2 and 14.9
    // mm.  The circle's 5.072 mm is stable to a thousandth over the same
    // change.  Three times is the claim; a fourth significant figure on the
    // square would be a claim about frame alignment.
    ASSERT_TRUE(q.max_error > 2.0 * c.max_error);
}

}  // namespace

int main() {
    test_every_offered_setting_keeps_the_ball();
    test_feedforward_is_worth_having();
    test_a_corner_is_harder_than_a_curve();
    std::printf("test_trajectory: all passed\n");
    return 0;
}
