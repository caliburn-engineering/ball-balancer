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
#include "sim_step.h"
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

// The opening, run through `stepSim` — the same step `PlateView::step` drives,
// and the reason this file is evidence about the application rather than about
// itself.  The setpoint driven before the loop reads it, the airborne ball
// shown where it will land, the servo rate taken at the instant the pose is
// for: all of that is inside the step now, in one place, and this harness
// cannot disagree with the demo about any of it.  See #30.
//
// What is left here is the measuring: a sample per frame of where the ball was
// when the frame opened, against the setpoint that frame was chasing.
DemoRun runDemo(const ModelEntry& e, const Eigen::MatrixXd& K,
                const SetpointPath& path, double duration) {
    const SimPlate plate = cascadePlate(e.params);

    SimInput in;
    in.design = cascadeDesign(e.params);
    in.design.K = K;
    in.path = path;

    SimState s = simStart(plate, in.design.home_leg_rad, attractStart(path));

    const double dt = in.dt;
    const int steps = static_cast<int>(duration / dt);

    DemoRun r;
    double t = 0.0;
    for (int k = 0; k < steps; ++k) {
        // Where the ball and the legs were when the frame opened — the sample
        // this trace is made of, taken before the step moves either.
        const Eigen::Matrix<double, 6, 1> open =
            plateFrame(s.ball, s.motion, kBallRadius);
        const double span = *std::max_element(s.alpha_rad.begin(), s.alpha_rad.end()) -
                            *std::min_element(s.alpha_rad.begin(), s.alpha_rad.end());

        const SimReport frame = stepSim(plate, in, s);
        const Eigen::Vector2d sp = frame.setpoint;

        r.trace.push_back({t, open(0), open(1), std::hypot(open(0), open(1)),
                           std::hypot(open(0) - sp(0), open(1) - sp(1)), span});
        t += dt;

        if (frame.airborne) ++r.airborne_frames;
        r.peak_height = std::max(r.peak_height, frame.ball_plate(2) - kBallRadius);
        if (frame.left_plate) {
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

    /// The least the plate ever pressed the ball, in m/s^2 — `N/m`, the very
    /// quantity separation is decided by, taken from the step rather than
    /// recomputed out here.  A separation count alone cannot say whether a
    /// tuning is comfortably holding the ball or missing by a hair's breadth,
    /// and those are different claims about the same demo.
    ///
    /// Only meaningful over a sweep that never separated: `SimReport` reports
    /// zero for a ball already in flight, which would drag this to zero and
    /// say nothing.  The one test that reads it asserts `separated == 0` first.
    double min_normal = 1e9;
};

// One shove, from the state the loop actually settles into — the ball parked a
// few millimetres off centre in the friction dead band, the legs slightly off
// home.  Shoving a pristine plate at exact centre would be an easier test than
// the thing it stands in for: #22's real failure had the ball at +3.4, -3.3 mm
// with the legs already displaced.
Sweep sweepOneDirection(const ModelEntry& e, const Eigen::MatrixXd& K,
                        const Disturbance& s, double theta) {
    const SimPlate plate = cascadePlate(e.params);

    SimInput in;
    in.design = cascadeDesign(e.params);
    in.design.K = K;
    // A held setpoint at the centre, which is `SimInput`'s default: the ticket
    // asks whether the loop brings the ball home, not whether it can track.

    SimState st = simStart(plate, in.design.home_leg_rad,
                           Eigen::Vector4d(s.start_x, s.start_y, 0.0, 0.0));

    const double dt = in.dt;
    Sweep sw;
    double peak = 0.0;
    Eigen::Matrix<double, 6, 1> end = Eigen::Matrix<double, 6, 1>::Zero();
    // Settle to the centre first, then shove from wherever that left everything.
    const int settle = static_cast<int>(s.settle_s / dt);
    const int total = settle + static_cast<int>(12.0 / dt);
    for (int k = 0; k < total; ++k) {
        // The shove: a step in the ball's own velocity, applied between two
        // frames.  The only thing this harness does to the simulation that the
        // Nudge buttons do not.
        if (k == settle && !st.ball.airborne) {
            st.ball.rolling(2) += s.speed * std::cos(theta);
            st.ball.rolling(3) += s.speed * std::sin(theta);
        }

        const SimReport frame = stepSim(plate, in, st);

        // The plate always has an assembly now.  A stronger statement than
        // "the ball stayed on", and the one #22's fix actually makes.
        ASSERT_TRUE(frame.fk.converged);
        ASSERT_TRUE(!frame.left_plate);

        if (frame.airborne) sw.separated = true;
        end = frame.ball_plate;
        sw.min_normal = std::min(sw.min_normal, frame.normal_accel);
        if (k >= settle) peak = std::max(peak, std::hypot(end(0), end(1)));
        sw.peak_height = std::max(sw.peak_height, end(2) - kBallRadius);
    }
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
    // 30.7 mm measured over these 720 directions, against the 280 mm the plate
    // has.  This bound was raised to 90 mm when the ball first learned to hop,
    // on the reasoning that a hopping ball carries its speed without rolling
    // friction and so lands further out — true in general, and no longer what
    // happens here: the shipped tuning does not separate the ball at this
    // disturbance at all (see the test below).  The bound is left where it is
    // rather than tightened onto the measurement, because what it is for is
    // catching a recovery that swings toward the rim.
    ASSERT_TRUE(worst_peak < 0.09);
}

// #23 asks for the shipped tuning's behaviour to be STATED: does it hop, and
// how often.  Stating it in CONTEXT.md makes it prose that can rot; stating it
// here makes it a fact that fails when it stops being true.
//
// **The answer has moved twice, and both times because the measurement was
// being taken against something other than the plate the demo runs.**  It said
// "30 of 72 directions" and it is now none of them:
//
//   - Those 30 were measured while the plate's accelerations were a finite
//     difference of two frames.  The leg command is a zero-order hold, so
//     differencing the rates it drives answers a command step with a delta
//     function — most of those hops were the estimator and not the plate.
//     Making the accelerations analytic (#23) took the count at this
//     disturbance to zero on its own.
//   - What survived was measured by a harness that took the servo rate at the
//     frame's START while pairing it with the pose at the frame's END, which
//     runs the plate 40 per cent fast.  Driving the application's own step
//     removed the rest (#30): 26 of 72 directions separated at 0.30 m/s under
//     the old harness, none do now.
//
// So, measured against the plate the application actually steps: **the shipped
// tuning does not throw the ball** at the disturbance the demo's own Nudge
// buttons offer.  Not marginally — the worst frame of 72 directions still has
// 1.62 m/s^2 of normal force in hand, a sixth of a g.
//
// What this deliberately does NOT do is map where the tuning DOES start to let
// go.  That is a sweep over speed as well as direction, it is
// [#31](https://github.com/caliburn-engineering/caliburn/issues/31)'s, and it
// is one of the harnesses #30 exists to be written against.  The one-line
// answer, unpinned until then: above about 0.30 m/s, in narrow slivers of
// direction, by fractions of a millimetre.
void test_the_shipped_tuning_holds_the_ball_through_its_own_disturbance() {
    const auto models = getBuiltinModels();
    const auto& e = cascadeModel(models);
    const Eigen::MatrixXd K = defaultGain(e);

    double margin = 1e9;
    int separated = 0;
    for (int i = 0; i < 72; ++i) {
        const Sweep sw = sweepOneDirection(e, K, Disturbance{},
                                           i * 2.0 * M_PI / 72.0);
        if (sw.separated) ++separated;
        margin = std::min(margin, sw.min_normal);
    }
    ASSERT_EQ(separated, 0);
    // And with room to spare rather than by a whisker: 1.62 measured.  The bar
    // is 1.0, which still fails if the plate ever comes within a tenth of a g
    // of letting go at the disturbance the demo hands out.  Asserted only
    // because `separated` is zero — see `Sweep::min_normal`.
    ASSERT_TRUE(margin > 1.0);
}

// ---------------------------------------------------------------------------
// A shove landing on a plate that is already working
// ---------------------------------------------------------------------------
//
// Every sweep above sets the ball down at rest, lets it settle to the centre of
// a level plate against a HELD setpoint, and shoves it there.  **The demo never
// looks like that.**  It opens tracking a circle and stays tracking it, so when
// a visitor reaches for Nudge the legs are already displaced, already moving,
// and the setpoint is somewhere else on the lap.  The shove lands on a plate
// that is already spending its workspace.
//
// That gap is why the interface could offer a shove the loop could not take
// while every test passed.  Measured over 36 points of the lap and 24
// directions, adding the shove to a ball tracking the opening circle:
//
//     speed   Nominal   Aggressive   Detuned      (failures out of 864)
//     0.25       0           0          0
//     0.30       0           0          0
//     0.35       6           5          0
//     0.40       9           5          1
//     0.50      37         123          0
//
// From rest all three survive 0.5 in every direction, which is what
// `kMaxNudgeSpeed` used to be set to and what
// `test_aggressive_is_no_more_fragile_than_the_shipped_tuning` still checks.
// While tracking, 0.5 loses the ball at every point of the lap once the two
// Nudge buttons are pressed together.  So the bound came down to 0.30, and the
// slider to `kMaxNudgePerAxis` so the pair composes to exactly that.
//
// Detuned's row is the reminder that these are slivers rather than a threshold:
// it fails at 0.40 and 0.45 and is clean again at 0.50.  The bound is set below
// the first failure of any tuning, with margin, exactly as `kMaxSetpointSpeed`
// is.
struct Tracked {
    bool lost = false;
    double end_error = 0.0;   ///< [m] ball to setpoint, twenty seconds later
};

/// One shove of `speed` in direction `theta`, delivered `phase0` of the way
/// round the lap while the opening circle is being tracked.
Tracked shoveWhileTracking(const ModelEntry& e, const Eigen::MatrixXd& K,
                           double phase0, double speed, double theta) {
    const SimPlate plate = cascadePlate(e.params);
    SimInput in;
    in.design = cascadeDesign(e.params);
    in.design.K = K;
    in.path = openingPath();

    SimState s = simStart(plate, in.design.home_leg_rad, attractStart(in.path));

    const double dt = in.dt;
    const int at = static_cast<int>(phase0 * in.path.period_s / dt);
    const int total = at + static_cast<int>(20.0 / dt);

    Tracked r;
    for (int k = 0; k < total; ++k) {
        // Exactly what the buttons do: added to the ball's own velocity, not
        // assigned over it.  The ball is already running at 75 mm/s, and a
        // shove is a shove rather than a teleport of the state.
        if (k == at && !s.ball.airborne) {
            s.ball.rolling(2) += speed * std::cos(theta);
            s.ball.rolling(3) += speed * std::sin(theta);
        }
        const SimReport f = stepSim(plate, in, s);
        if (f.left_plate) { r.lost = true; return r; }
        if (k + 1 == total)
            r.end_error = std::hypot(f.ball_plate(0) - f.setpoint(0),
                                     f.ball_plate(1) - f.setpoint(1));
    }
    return r;
}

// The hardest shove the interface can compose, at every point of the lap and
// from every direction: the ball stays on the plate and goes back to tracking.
//
// A coarser grid than the 36 x 24 the table above was measured on, because this
// runs on every build; the bound it checks is the one that measurement chose.
void test_a_shove_while_tracking_is_rejected_from_every_direction() {
    const auto models = getBuiltinModels();
    const auto& e = cascadeModel(models);
    const Eigen::MatrixXd K = defaultGain(e);

    double worst_error = 0.0;
    for (int i = 0; i < 12; ++i) {
        for (int j = 0; j < 12; ++j) {
            const Tracked t = shoveWhileTracking(
                e, K, double(i) / 12.0, kMaxNudgeSpeed, j * M_PI / 6.0);
            ASSERT_TRUE(!t.lost);
            worst_error = std::max(worst_error, t.end_error);
        }
    }
    // And back ON the path, not merely still on the plate.  The demo's own
    // tracking error is 3.7 mm; twenty seconds after the worst shove the
    // interface can deliver it is back inside 15.
    ASSERT_TRUE(worst_error < 0.015);
}

// The buttons compose, and the bound has to survive that.
//
// "Nudge +x" and "Nudge +y" each add the slider's value to one axis, so the
// worst the pair can do is `sqrt(2)` times the slider at 45 degrees.  The
// slider's top used to be `kMaxNudgeSpeed` itself, which made the constant's
// own description — the largest disturbance the interface can hand the loop —
// false by a factor of 1.41, and handed it 0.707 m/s.  Arithmetic, so it is
// checked as arithmetic: raise the slider's ceiling without re-measuring the
// envelope and this fails.
void test_the_nudge_buttons_cannot_compose_past_the_bound() {
    const double both = std::hypot(kMaxNudgePerAxis, kMaxNudgePerAxis);
    ASSERT_NEAR(both, kMaxNudgeSpeed, 1e-12);
    ASSERT_TRUE(kMaxNudgePerAxis < kMaxNudgeSpeed);
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
//
// > **This is currently failing, and it is not #30's to fix.**
// > Direction 33 of 180 — 66.0 degrees — rolls the ball to 293.7 mm against a
// > rim at 280 mm and loses it.  Not a hop: `separated` is false the whole way,
// > so the aggressive gain simply flings it wider than the plate.  The same
// > disease under the harness's old servo-rate instant showed up at a different
// > grid point instead (the shipped tuning, 1 of 720 directions, at 162.5
// > degrees), which is what `test_attract_mode` was already failing on before
// > this change; the last tree where every direction held was `3303ff4`, one
// > commit before the plate's accelerations became analytic.  So this is #23's
// > to answer — either the analytic normal force needs a further look, or the
// > Aggressive preset's `Q` of 150 no longer clears the bar it was chosen for.
void test_the_aggressive_preset_recovers_from_every_direction() {
    const auto models = getBuiltinModels();
    const auto& e = cascadeModel(models);

    const Eigen::MatrixXd K = gainForPreset(e, presetNamed("Aggressive"));

    for (int i = 0; i < 180; ++i)
        sweepOneDirection(e, K, Disturbance{}, i * M_PI / 90.0);
}

// And the fact that only exists now the ball can leave: an over-aggressive
// tuning does not merely overshoot, it takes the ball OFF THE SURFACE and then
// off the plate.  Q at 2000 slams the legs hard enough that the plate is
// moving out from under the ball rather than tilting under it.
//
// Measured against the application's own step, at the demo's own 0.26 m/s
// disturbance and over 90 directions: the ball leaves the surface in 4 of them
// and leaves the plate entirely in 1.  Before #23 the same tuning looked
// merely fast — 0.35 s to settle — because a ball glued to the plate cannot be
// thrown off it.
//
// **The size of the hop is nothing like what this test used to claim.**  It
// said "one launch reached 689 mm of altitude"; the ball now rises 1.5 mm.
// Both of the reasons are in
// `test_the_shipped_tuning_holds_the_ball_through_its_own_disturbance` above —
// a differenced acceleration answering a stepping command with a delta
// function (#23), and a harness taking the servo rate at the wrong instant
// (#30).  A metre of altitude off a 0.26 m/s nudge was never physics.
//
// What survives, and is the point of the test, is that losing the ball is a
// thing an over-aggressive gain does and a sensible one does not.  That is
// still the honest ceiling on #19's aggressive preset, and still a far better
// demonstration of over-aggressive control than a settling time.
void test_an_over_aggressive_tuning_throws_the_ball_off() {
    const auto models = getBuiltinModels();
    const auto& e = cascadeModel(models);

    Eigen::VectorXd q(7);
    q << 1, 1, 1, 2000, 2000, 2, 2;
    const Eigen::MatrixXd K = gainFor(e, q, defaultLqrInputWeights(3));

    const SimPlate plate = cascadePlate(e.params);
    SimInput in;
    in.design = cascadeDesign(e.params);
    in.design.K = K;

    const Disturbance s;
    const double dt = in.dt;
    const int settle = static_cast<int>(s.settle_s / dt);

    int lost = 0;
    bool separated = false;
    double peak_hop = 0.0;
    for (int i = 0; i < 90; ++i) {
        const double theta = i * M_PI / 45.0;
        SimState st = simStart(plate, in.design.home_leg_rad,
                               Eigen::Vector4d(s.start_x, s.start_y, 0.0, 0.0));

        for (int k = 0; k < settle + static_cast<int>(6.0 / dt); ++k) {
            if (k == settle && !st.ball.airborne) {
                st.ball.rolling(2) += s.speed * std::cos(theta);
                st.ball.rolling(3) += s.speed * std::sin(theta);
            }
            const SimReport frame = stepSim(plate, in, st);
            if (frame.airborne) separated = true;
            peak_hop = std::max(peak_hop, frame.ball_plate(2) - kBallRadius);
            if (frame.left_plate) { ++lost; break; }
        }
    }
    ASSERT_TRUE(lost > 0);            // it really does lose the ball
    ASSERT_TRUE(separated);           // and really does take it off the surface
    // Millimetres, not metres: 1.5 mm measured over these 90 directions.  Both
    // bounds are assertions — too small and the plate has stopped letting go at
    // all, too large and something is differencing a step again.
    ASSERT_TRUE(peak_hop > 0.0005);
    ASSERT_TRUE(peak_hop < 0.020);
}

}  // namespace

int main() {
    test_the_demo_opens_on_the_path_and_moving_with_it();
    test_the_opening_path_is_a_gentle_circle();
    test_the_opening_is_moving_within_a_second();
    test_the_opening_does_not_slam_the_legs();
    test_the_demo_tracks_its_circle();
    test_the_kick_is_rejected_from_every_direction();
    test_the_shipped_tuning_holds_the_ball_through_its_own_disturbance();
    test_a_shove_while_tracking_is_rejected_from_every_direction();
    test_the_nudge_buttons_cannot_compose_past_the_bound();
    test_the_aggressive_preset_recovers_from_every_direction();
    test_an_over_aggressive_tuning_throws_the_ball_off();
    test_ten_minutes_unattended_never_loses_the_ball();
    std::printf("test_attract_mode: all passed\n");
    return 0;
}
