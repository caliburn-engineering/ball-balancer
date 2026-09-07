// tests/test_attract_mode.cpp
//
// The demo running itself for a visitor who has not arrived yet, and every one
// of its acceptance criteria is phrased as something a person would see:
// "motion is apparent", "the controller is already active at load", "the ball
// never leaves the plate".  None of that can be checked by looking at a canvas
// once and calling it done — the tuning, the plate geometry and the rolling
// model all move underneath it.
//
// So the opening is a pure function of the path and of elapsed time, and this
// file runs it against the same nonlinear plate the browser does.  "Visible"
// becomes millimetres and "the visitor never sees a broken demo" becomes: the
// ball does not leave the plate in ten minutes of unattended running.
//
// **The demo used to kick the ball and does not any more.**  See
// `attract_mode.h` for why — briefly, a kick reads as a fault rather than as a
// disturbance when there is nothing still to read it against, and the opening
// displacement in particular slammed the legs 20.7 degrees apart in three
// frames.  Two consequences for this file:
//
//   - The opening tests are now about a ball tracing a circle from the first
//     frame, and one of them pins the leg swing that used to be the glitch.
//   - The disturbance sweeps stayed.  What they test is a property of the
//     plant and the gain — that a 0.26 m/s shove is rejected from every
//     direction without losing the ball — which is still true, still bounds
//     what #19 may offer, and is no less worth pinning because the demo no
//     longer performs it.  The kick moved from production into `Disturbance`
//     below, where it is honestly a test fixture.
#include "attract_mode.h"

#include "auto_balance.h"
#include "ball_contact.h"
#include "ball_sim.h"
#include "cascade_fixture.h"
#include "setpoint_path.h"
#include "table_kinematics.h"
#include "test_helpers.h"

#include <algorithm>
#include <cstdio>
#include <cmath>
#include <vector>

using namespace caliburn;

namespace {

constexpr double kBallRadius = kFixtureBallRadius;
constexpr double kDeg = M_PI / 180.0;

// "Home again", for a plant whose rolling friction makes the centre a REGION
// rather than a point.  A state feedback has no integral term, so it parks
// wherever the tilt it is asking for falls inside the dead band and then
// creeps.  10 mm on a 300 mm radius is a handful of pixels, which is the sense
// in which the ticket asks for a recovery to be visible.
constexpr double kHome = 0.010;  // [m]

// The disturbance the sweeps use, and the reason it is here rather than in
// `attract_mode`.
//
// These numbers were attract mode's own: the ball placed 72 mm out, left 2.5 s
// to settle, then shoved at 0.26 m/s.  The demo no longer does any of that.
// But "this gain rejects a 0.26 m/s shove from every direction" is a claim
// about the plant and the tuning, and #22 is the story of what happens when
// nobody checks it — so the numbers live on as a fixture.
//
// 0.26 m/s specifically, and not more: above roughly 0.29 there are narrow
// slivers of direction, a few in every thousand, where the gain saturates
// against the workspace and cannot bring the ball back.  0.28 loses 2
// directions in 5760, 0.29 loses 6, 0.26 loses none at that resolution.  The
// slivers thin as the kick shrinks rather than vanishing at a threshold, and
// nobody has proved there is not a narrower one still.
struct Disturbance {
    double start_x = 0.06;    ///< [m] where the ball is placed to settle from
    double start_y = -0.04;   ///< [m]
    double settle_s = 2.5;    ///< [s] before the shove
    double speed = 0.26;      ///< [m/s]
};

// ---------------------------------------------------------------------------
// The opening
// ---------------------------------------------------------------------------

// The demo opens ON its path and MOVING ALONG it.  Both halves, because the
// reference state carries the path's velocity: a ball placed correctly but at
// rest still hands the loop a 75 mm/s error to answer on the first frame, and
// answering it is the slam this replaced.
void test_the_demo_opens_on_the_path_and_moving_with_it() {
    const SetpointPath p = openingPath();
    const Eigen::Vector4d b = attractStart(p);

    const Eigen::Vector2d p0 = pathPoint(p, 0.0);
    const Eigen::Vector2d v0 = pathVelocity(p, 0.0);
    ASSERT_NEAR(b(0), p0(0), 1e-15);
    ASSERT_NEAR(b(1), p0(1), 1e-15);
    ASSERT_NEAR(b(2), v0(0), 1e-15);
    ASSERT_NEAR(b(3), v0(1), 1e-15);

    // And it is genuinely moving, not merely "on a path" that happens to be
    // stationary at t = 0.
    ASSERT_TRUE(std::hypot(b(2), b(3)) > 0.05);
}

// The opening path is a circle, and a gentle one: the corners are the better
// demonstration but they are an argument the visitor should choose to hear.
void test_the_opening_path_is_a_gentle_circle() {
    const SetpointPath p = openingPath();
    ASSERT_TRUE(p.shape == PathShape::Circle);

    // Comfortably inside both bounds rather than up against either — the
    // opening's job is to look effortless.
    const double speed = pathLength(p) / p.period_s;
    ASSERT_TRUE(speed < 0.4 * kMaxSetpointSpeed);
    ASSERT_TRUE(p.radius_m < 0.75 * kMaxPathRadius);
    ASSERT_TRUE(p.period_s >= minPeriod(p));

    // And well inside the plate, with room for the ball to trail and overshoot.
    ASSERT_TRUE(p.radius_m < 0.5 * 0.30);
}

// ---------------------------------------------------------------------------
// The demo, run against the nonlinear plate
// ---------------------------------------------------------------------------

struct Sample {
    double t;
    double x, y;      ///< [m] ball, plate frame
    double radius;    ///< [m] from the centre
    double error;     ///< [m] from the setpoint the ball is chasing
    double leg_span;  ///< [rad] widest minus narrowest leg, this frame
};

struct DemoRun {
    std::vector<Sample> trace;
    bool left_plate = false;
    int airborne_frames = 0;   ///< the ball can leave upward now, see #23
    double peak_height = 0.0;  ///< [m] above the plate surface
};

// PlateView::step's calls, in PlateView::step's order: drive the setpoint from
// the path, legCommand -> stepServosOnPlate -> solve_pose -> stepBallContact.
// The setpoint is driven BEFORE the loop reads it, so the gain sees this
// frame's target rather than last frame's — the same order `plate_view` runs,
// which is what makes this test evidence about that code.
DemoRun runDemo(const ModelEntry& e, const Eigen::MatrixXd& K,
                const SetpointPath& path, double duration) {
    const TableParams tp = cascadeMechanism(e.params);
    const TableKinematics tk(tp);

    AutoBalanceDesign d;
    d.K = K;
    d.home_leg_rad = cascadeHomeLegAngle(e.params);
    d.servo_tau = cascadeServoTau(e.params);
    d.alpha_min_rad = tp.alpha_min;
    d.alpha_max_rad = tp.alpha_max;

    const RollingBallDynamics dynamics(
        kPlateBall, PlateParams{tp.R_table, cascadeGravity(e.params)});

    std::array<double, 3> alpha = {d.home_leg_rad, d.home_leg_rad, d.home_leg_rad};
    BallState ball;
    ball.rolling = attractStart(path);
    TablePose pose = tk.home_pose(d.home_leg_rad);
    PlateMotion pm = plateMotion(tk, pose, alpha, {0.0, 0.0, 0.0}), pm_prev = pm;

    const double dt = 1.0 / 60.0;
    const int steps = static_cast<int>(duration / dt);
    const double g = cascadeGravity(e.params);

    DemoRun r;
    double t = 0.0;
    // The phase the application drives, through the function the application
    // drives it with.  Deriving it from `t` is the shape of the bug #24 fixed,
    // and a harness that keeps the old shape stops being evidence about this.
    double phase = 0.0;
    for (int k = 0; k < steps; ++k) {
        const PathStep step = stepPath(path, phase, dt);
        const Eigen::Vector2d sp = step.point;

        const Eigen::Matrix<double, 6, 1> bp = plateFrame(ball, pm, kBallRadius);
        Eigen::Vector4d seen(bp(0), bp(1), bp(3), bp(4));
        if (ball.airborne) {
            const Eigen::Vector2d land = predictedLanding(bp, kBallRadius, g);
            seen << land(0), land(1), 0.0, 0.0;
        }
        // The reference the loop is given: the setpoint, and — while the ball
        // is on the plate — the speed the setpoint is travelling at.  Zero in
        // the air, where the plate cannot touch the ball and `seen` is already
        // the predicted landing point.  Exactly what `plate_view` builds.
        BallReference ref;
        ref.position = sp;
        if (!ball.airborne) ref.velocity = step.velocity;

        const double span = *std::max_element(alpha.begin(), alpha.end()) -
                            *std::min_element(alpha.begin(), alpha.end());
        r.trace.push_back({t, bp(0), bp(1), std::hypot(bp(0), bp(1)),
                           std::hypot(bp(0) - sp(0), bp(1) - sp(1)), span});

        const LegCommand c = legCommand(tk, d, alpha, seen, ref);
        t += dt;
        phase = step.next_phase;
        std::array<double, 3> adot{};
        for (int i = 0; i < 3; ++i)
            adot[i] = (c.alpha_rad[i] - alpha[i]) / d.servo_tau;
        alpha = stepServosOnPlate(tk, alpha, c.alpha_rad, d.servo_tau, dt);

        const FKResult fk = tk.solve_pose(alpha, pose);
        if (fk.converged) pose = fk.pose;
        pm_prev = pm;
        pm = plateMotion(tk, pose, alpha, adot,
                         servoAccel(adot, d.servo_tau), &pm_prev, dt);

        ball = stepBallContact(dynamics, ball, pm, pm_prev, pose,
                               kBallRadius, g, dt);
        if (ball.airborne) ++r.airborne_frames;

        const Eigen::Matrix<double, 6, 1> now = plateFrame(ball, pm, kBallRadius);
        r.peak_height = std::max(r.peak_height, now(2) - kBallRadius);
        if (!ballOnPlate(Eigen::Vector4d(now(0), now(1), now(3), now(4)),
                         tp.R_table, kBallRadius)) {
            r.left_plate = true;
            return r;
        }
    }
    return r;
}

double peakErrorAfter(const DemoRun& r, double t0) {
    double peak = 0.0;
    for (const Sample& s : r.trace)
        if (s.t >= t0) peak = std::max(peak, s.error);
    return peak;
}

double meanErrorAfter(const DemoRun& r, double t0) {
    double sum = 0.0;
    int n = 0;
    for (const Sample& s : r.trace)
        if (s.t >= t0) { sum += s.error; ++n; }
    return n ? sum / n : 0.0;
}

double peakSpanBefore(const DemoRun& r, double t1) {
    double peak = 0.0;
    for (const Sample& s : r.trace) {
        if (s.t > t1) break;
        peak = std::max(peak, s.leg_span);
    }
    return peak;
}

// "Motion is apparent within the first few seconds."  As a number: the ball is
// tracing a 120 mm circle at 75 mm/s, so in one second it sweeps about a tenth
// of a lap — 36 degrees, some 75 mm of travel.  A still image cannot do that,
// and neither can a demo that waits for a click.
//
// Measured as angle swept about the centre, not as distance from it: a circle
// keeps the radius constant on purpose, so the old test — which watched the
// radius fall as the ball came home from its opening displacement — is
// measuring a motion that no longer exists.
void test_the_opening_is_moving_within_a_second() {
    const auto models = getBuiltinModels();
    const auto& e = cascadeModel(models);

    const DemoRun r = runDemo(e, defaultGain(e), openingPath(), 2.0);
    ASSERT_TRUE(!r.left_plate);

    const Sample& first = r.trace.front();
    const Sample* at_one = nullptr;
    for (const Sample& s : r.trace)
        if (s.t >= 1.0) { at_one = &s; break; }
    ASSERT_TRUE(at_one != nullptr);

    double swept = std::atan2(at_one->y, at_one->x) - std::atan2(first.y, first.x);
    while (swept < 0.0) swept += 2.0 * M_PI;
    // A tenth of a lap, less a margin for the ball trailing the setpoint.
    ASSERT_TRUE(swept > 25.0 * kDeg);

    // And it is out on the circle throughout, not drifting in to the centre.
    ASSERT_TRUE(first.radius > 0.05);
    ASSERT_TRUE(at_one->radius > 0.05);
}

// The glitch, as a regression test.
//
// This is the one the ticket was actually about.  Opening with the ball 72 mm
// off centre handed the state feedback a step input at t = 0, with the legs
// exactly at home and no servo history to smear it: measured, the legs swung
// **20.7 degrees apart within three frames** and the table dropped 4.5 mm
// before settling inside 250 ms.  Correct, and it read as the mechanism
// glitching.
//
// Starting on the path makes both errors zero at t = 0, and the same
// measurement gives 3.0 degrees.  The bar is set at 8, which is comfortably
// under the old behaviour and comfortably over the new one: what must never
// come back is a demo whose first visible act is a spasm.
void test_the_opening_does_not_slam_the_legs() {
    const auto models = getBuiltinModels();
    const auto& e = cascadeModel(models);

    const DemoRun r = runDemo(e, defaultGain(e), openingPath(), 3.0);
    ASSERT_TRUE(!r.left_plate);
    ASSERT_TRUE(peakSpanBefore(r, 1.0) < 8.0 * kDeg);

    // And the ball never leaves the plate on the opening frames — a hop in the
    // first half second would read as a glitch just as surely.
    ASSERT_EQ(r.airborne_frames, 0);
}

// The demo tracks its circle, and this is the number that says how well.
//
// 3.7 mm of mean error measured, on a 120 mm circle — the ball follows the
// path rather than sitting inside it, which is what the velocity feedforward
// bought (#24: 3.66 mm with it, 35.52 mm without).  The bar is set at 15 mm:
// loose enough not to be a transcript of one LQR solve, tight enough that
// losing the feedforward fails it by more than a factor of two.
void test_the_demo_tracks_its_circle() {
    const auto models = getBuiltinModels();
    const auto& e = cascadeModel(models);

    const DemoRun r = runDemo(e, defaultGain(e), openingPath(), 30.0);
    ASSERT_TRUE(!r.left_plate);

    // After the first lap, so this measures tracking rather than any opening
    // transient — though the whole point of the opening is that there is none.
    ASSERT_TRUE(meanErrorAfter(r, 10.0) < 0.015);
    ASSERT_TRUE(peakErrorAfter(r, 10.0) < 0.030);
}

// Ten minutes unattended.  The failure this pins is the one that shipped: the
// ball left at 59.3 s and then left on every kick after it, for as long as the
// page stayed open, because the plate had folded onto an assembly it could
// never climb out of.
//
// Sixty laps of the circle now rather than 150 kicks, which is a different and
// arguably harder question: the old run spent most of its time with the plate
// nearly level between disturbances, and this one never stops steering.
void test_ten_minutes_unattended_never_loses_the_ball() {
    const auto models = getBuiltinModels();
    const auto& e = cascadeModel(models);

    const DemoRun r = runDemo(e, defaultGain(e), openingPath(), 600.0);
    ASSERT_TRUE(!r.left_plate);

    // And never even came close: the ball leaves at R_table - r_ball, and
    // tracking a 120 mm circle should never take it past about 140.
    const TableParams tp = cascadeMechanism(e.params);
    double peak_radius = 0.0;
    for (const Sample& s : r.trace) peak_radius = std::max(peak_radius, s.radius);
    ASSERT_TRUE(peak_radius < 0.5 * (tp.R_table - kBallRadius));

    // The error does not creep over ten minutes.  A slow drift would be
    // invisible in a thirty-second test and obvious to a visitor watching a
    // kiosk, which is exactly the failure this length of run is for.
    ASSERT_TRUE(meanErrorAfter(r, 570.0) < 0.015);

    // Tracking a gentle circle does not throw the ball.  The shipped tuning
    // CAN hop — see the sweep below — but only when something shoves it, and
    // an unattended demo now has nothing that does.
    ASSERT_EQ(r.airborne_frames, 0);
}

// ---------------------------------------------------------------------------
// Disturbance rejection
// ---------------------------------------------------------------------------
//
// The demo no longer performs this, and the loop had better still do it: the
// visitor has Nudge buttons, #19 will offer tunings that change the answer,
// and #22 was a permanent failure that a sweep like this would have caught on
// the day it landed.

/// What one shove did: how far out the ball got, and whether it left the plate
/// on the way.  Separation is reported rather than inferred from the reach,
/// because a hop and a wide roll look identical in a radius.
struct Sweep {
    double peak_radius = 0.0;   ///< [m] from centre
    bool separated = false;
    double peak_height = 0.0;   ///< [m] above the surface
};

// One shove, from the state the loop actually settles into — the ball parked a
// few millimetres off centre in the friction dead band, the legs slightly off
// home.  Shoving a pristine plate at exact centre would be an easier test than
// the thing it stands in for: #22's real failure had the ball at +3.4, -3.3 mm
// with the legs already displaced.
Sweep sweepOneDirection(const ModelEntry& e, const Eigen::MatrixXd& K,
                        const Disturbance& s, double theta) {
    const TableParams tp = cascadeMechanism(e.params);
    const TableKinematics tk(tp);
    AutoBalanceDesign d;
    d.K = K;
    d.home_leg_rad = cascadeHomeLegAngle(e.params);
    d.servo_tau = cascadeServoTau(e.params);
    d.alpha_min_rad = tp.alpha_min;
    d.alpha_max_rad = tp.alpha_max;
    const RollingBallDynamics dynamics(
        kPlateBall, PlateParams{tp.R_table, cascadeGravity(e.params)});

    std::array<double, 3> alpha = {d.home_leg_rad, d.home_leg_rad, d.home_leg_rad};
    TablePose pose = tk.home_pose(d.home_leg_rad);
    BallState ball;
    ball.rolling << s.start_x, s.start_y, 0.0, 0.0;
    PlateMotion pm = plateMotion(tk, pose, alpha, {0.0, 0.0, 0.0}), pm_prev = pm;
    const double g = cascadeGravity(e.params);

    const double dt = 1.0 / 60.0;
    Sweep sw;
    double peak = 0.0;
    // Settle to the centre first, then shove from wherever that left everything.
    const int settle = static_cast<int>(s.settle_s / dt);
    const int total = settle + static_cast<int>(12.0 / dt);
    for (int k = 0; k < total; ++k) {
        if (k == settle && !ball.airborne) {
            ball.rolling(2) += s.speed * std::cos(theta);
            ball.rolling(3) += s.speed * std::sin(theta);
        }
        const Eigen::Matrix<double, 6, 1> bp = plateFrame(ball, pm, kBallRadius);
        Eigen::Vector4d seen(bp(0), bp(1), bp(3), bp(4));
        if (ball.airborne) {
            const Eigen::Vector2d land = predictedLanding(bp, kBallRadius, g);
            seen << land(0), land(1), 0.0, 0.0;
        }
        const LegCommand c = legCommand(tk, d, alpha, seen, {});
        std::array<double, 3> adot{};
        for (int i = 0; i < 3; ++i)
            adot[i] = (c.alpha_rad[i] - alpha[i]) / d.servo_tau;
        alpha = stepServosOnPlate(tk, alpha, c.alpha_rad, d.servo_tau, dt);
        const FKResult fk = tk.solve_pose(alpha, pose);
        // The plate always has an assembly now.  A stronger statement than
        // "the ball stayed on", and the one #22's fix actually makes.
        ASSERT_TRUE(fk.converged);
        pose = fk.pose;
        pm_prev = pm;
        pm = plateMotion(tk, pose, alpha, adot,
                         servoAccel(adot, d.servo_tau), &pm_prev, dt);
        ball = stepBallContact(dynamics, ball, pm, pm_prev, pose,
                               kBallRadius, g, dt);
        if (ball.airborne) sw.separated = true;
        const Eigen::Matrix<double, 6, 1> now = plateFrame(ball, pm, kBallRadius);
        if (k >= settle) peak = std::max(peak, std::hypot(now(0), now(1)));
        sw.peak_height = std::max(sw.peak_height, now(2) - kBallRadius);
        ASSERT_TRUE(ballOnPlate(Eigen::Vector4d(now(0), now(1), now(3), now(4)),
                                tp.R_table, kBallRadius));
    }
    const Eigen::Matrix<double, 6, 1> end = plateFrame(ball, pm, kBallRadius);
    ASSERT_TRUE(std::hypot(end(0), end(1)) < kHome);   // and came home
    sw.peak_radius = peak;
    return sw;
}

// Every direction, not the handful a schedule would reach.  A shove the loop
// cannot reject is one the demo's Nudge button must not offer, wherever it
// points.
void test_the_kick_is_rejected_from_every_direction() {
    const auto models = getBuiltinModels();
    const auto& e = cascadeModel(models);
    const Eigen::MatrixXd K = defaultGain(e);

    double worst_peak = 0.0;
    for (int i = 0; i < 720; ++i)
        worst_peak = std::max(worst_peak,
                              sweepOneDirection(e, K, Disturbance{},
                                                i * M_PI / 360.0).peak_radius);

    // Visible, and nowhere near the edge: 30 mm against the 280 mm the plate
    // has.  Measured; the assertion is loose around it on purpose.
    ASSERT_TRUE(worst_peak > 0.018);
    // Wider than it was before #23: a ball that hops carries its speed without
    // rolling friction, so it lands further out.  Still a quarter of what the
    // plate has.
    ASSERT_TRUE(worst_peak < 0.09);
}

// #23 asks for the shipped tuning's behaviour to be STATED: does it hop, and
// how often.  Stating it in CONTEXT.md makes it prose that can rot; stating it
// here makes it a fact that fails when it stops being true.
//
// 30 of 72 directions, measured — a bit under half, and the same 30 whether
// rolling resistance is scaled by the normal force or by g, because separation
// is decided by `normalAccel` and friction does not enter it.
//
// Note what this does and does not say now that the demo has stopped kicking.
// The shipped TUNING still hops when the ball is shoved, and a visitor pressing
// Nudge will see it.  The shipped DEMO no longer hops on its own, because
// nothing shoves the ball unless somebody asks — which the ten-minute test
// above pins from the other side.
void test_the_shipped_tuning_hops_in_a_stated_fraction_of_directions() {
    const auto models = getBuiltinModels();
    const auto& e = cascadeModel(models);
    const Eigen::MatrixXd K = defaultGain(e);

    int separated = 0;
    double worst_hop = 0.0;
    for (int i = 0; i < 72; ++i) {
        const Sweep sw = sweepOneDirection(e, K, Disturbance{},
                                           i * 2.0 * M_PI / 72.0);
        if (sw.separated) ++separated;
        worst_hop = std::max(worst_hop, sw.peak_height);
    }

    ASSERT_TRUE(separated > 15);
    ASSERT_TRUE(separated < 50);
    // And briefly: 6.6 mm measured, against a 40 mm ball.  The hop is meant to
    // be legible in the 3D view, not to dominate it.
    ASSERT_TRUE(worst_hop > 0.002);
    ASSERT_TRUE(worst_hop < 0.030);
}

// The Aggressive preset (#19), against the disturbance the demo itself
// offers.  It saturates the servos from every direction and still brings the
// ball home from every direction, which is what makes saturation a bounded
// cost rather than a failure — and what makes it honest to put the tuning
// behind a button a visitor is invited to press.
//
// `sweepOneDirection` asserts the ball stays on the plate and comes home, so
// the loop body is the whole test.  180 directions, not the handful a schedule
// would reach: a preset that fails in one direction is a preset that fails.
void test_the_aggressive_preset_recovers_from_every_direction() {
    const auto models = getBuiltinModels();
    const auto& e = cascadeModel(models);

    const Eigen::MatrixXd K = gainForPreset(e, presetNamed("Aggressive"));

    for (int i = 0; i < 180; ++i)
        sweepOneDirection(e, K, Disturbance{}, i * M_PI / 90.0);
}

// And the fact that only exists now the ball can leave: an over-aggressive
// tuning does not merely overshoot, it THROWS THE BALL OFF THE PLATE.
//
// Q at 2000 slams the legs hard enough that the plate is moving out from under
// the ball rather than tilting under it: measured, 54 of 720 directions end
// with the ball gone, and one launch reached 689 mm of altitude.  Before #23
// the same tuning looked merely fast — 0.35 s to settle — because a ball glued
// to the plate cannot be thrown off it.
//
// This is the honest ceiling on #19's aggressive preset, and a far better
// demonstration of over-aggressive control than a settling time.
void test_an_over_aggressive_tuning_throws_the_ball_off() {
    const auto models = getBuiltinModels();
    const auto& e = cascadeModel(models);

    Eigen::VectorXd q(7);
    q << 1, 1, 1, 2000, 2000, 2, 2;
    const Eigen::MatrixXd K = gainFor(e, q, defaultLqrInputWeights(3));

    const TableParams tp = cascadeMechanism(e.params);
    const TableKinematics tk(tp);
    AutoBalanceDesign d;
    d.K = K;
    d.home_leg_rad = cascadeHomeLegAngle(e.params);
    d.servo_tau = cascadeServoTau(e.params);
    d.alpha_min_rad = tp.alpha_min;
    d.alpha_max_rad = tp.alpha_max;
    const double g = cascadeGravity(e.params);
    const RollingBallDynamics dynamics(kPlateBall, PlateParams{tp.R_table, g});
    const Disturbance s;
    const double dt = 1.0 / 60.0;

    int lost = 0;
    double peak_hop = 0.0;
    for (int i = 0; i < 90; ++i) {
        const double theta = i * M_PI / 45.0;
        std::array<double, 3> alpha = {d.home_leg_rad, d.home_leg_rad, d.home_leg_rad};
        TablePose pose = tk.home_pose(d.home_leg_rad);
        BallState ball;
        ball.rolling << s.start_x, s.start_y, 0.0, 0.0;
        PlateMotion pm = plateMotion(tk, pose, alpha, {0.0, 0.0, 0.0}), pm_prev = pm;

        const int settle = static_cast<int>(s.settle_s / dt);
        for (int k = 0; k < settle + static_cast<int>(6.0 / dt); ++k) {
            if (k == settle && !ball.airborne) {
                ball.rolling(2) += s.speed * std::cos(theta);
                ball.rolling(3) += s.speed * std::sin(theta);
            }
            const Eigen::Matrix<double, 6, 1> bp = plateFrame(ball, pm, kBallRadius);
            Eigen::Vector4d seen(bp(0), bp(1), bp(3), bp(4));
            if (ball.airborne) {
                const Eigen::Vector2d land = predictedLanding(bp, kBallRadius, g);
                seen << land(0), land(1), 0.0, 0.0;
            }
            const LegCommand c = legCommand(tk, d, alpha, seen, {});
            std::array<double, 3> adot{};
            for (int j = 0; j < 3; ++j)
                adot[j] = (c.alpha_rad[j] - alpha[j]) / d.servo_tau;
            alpha = stepServosOnPlate(tk, alpha, c.alpha_rad, d.servo_tau, dt);
            const FKResult fk = tk.solve_pose(alpha, pose);
            if (fk.converged) pose = fk.pose;
            pm_prev = pm;
            pm = plateMotion(tk, pose, alpha, adot,
                         servoAccel(adot, d.servo_tau), &pm_prev, dt);
            ball = stepBallContact(dynamics, ball, pm, pm_prev, pose,
                                   kBallRadius, g, dt);
            const Eigen::Matrix<double, 6, 1> now = plateFrame(ball, pm, kBallRadius);
            peak_hop = std::max(peak_hop, now(2) - kBallRadius);
            if (!ballOnPlate(Eigen::Vector4d(now(0), now(1), now(3), now(4)),
                             tp.R_table, kBallRadius)) { ++lost; break; }
        }
    }
    ASSERT_TRUE(lost > 0);            // it really does lose the ball
    ASSERT_TRUE(peak_hop > 0.05);     // and really does launch it, metres not mm
}

}  // namespace

int main() {
    test_the_demo_opens_on_the_path_and_moving_with_it();
    test_the_opening_path_is_a_gentle_circle();
    test_the_opening_is_moving_within_a_second();
    test_the_opening_does_not_slam_the_legs();
    test_the_demo_tracks_its_circle();
    test_the_kick_is_rejected_from_every_direction();
    test_the_shipped_tuning_hops_in_a_stated_fraction_of_directions();
    test_the_aggressive_preset_recovers_from_every_direction();
    test_an_over_aggressive_tuning_throws_the_ball_off();
    test_ten_minutes_unattended_never_loses_the_ball();
    std::printf("test_attract_mode: all passed\n");
    return 0;
}
