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

    /// The worst `|cmd - alpha|` any leg was handed at a frame's start, in rad.
    ///
    /// This is the quantity the servo rate limit is a threshold on — it engages
    /// above `rate_max * tau`, 30 degrees against the shipped servo — so it is
    /// how a sweep says whether the limit was even in the room.  See
    /// `kServoRateMax` and #32.
    double peak_leg_error = 0.0;
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

        const std::array<double, 3> was = st.alpha_rad;
        const SimReport frame = stepSim(plate, in, st);
        for (int i = 0; i < 3; ++i)
            sw.peak_leg_error =
                std::max(sw.peak_leg_error, std::abs(frame.cmd_rad[i] - was[i]));

        // The plate always has an assembly now.  A stronger statement than
        // "the ball stayed on", and the one #22's fix actually makes.
        ASSERT_TRUE(frame.fk.converged);
        // And it is always the SAME assembly — the one the machine is built
        // in, which is the one `can_hold` vouched for when the command was
        // clipped.  #29 is a whole ticket about what happens when those two
        // part company: the plate came out of a near-singular frame on the
        // other assembly, tilted downhill toward the ball it was trying to
        // catch, and rolled it off.  Asked here rather than in one test,
        // because every sweep in this file is a chance for it.
        ASSERT_TRUE(onBuiltAssembly(plate.kinematics(), st.alpha_rad, st.pose));
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
    // 31.2 mm measured over these 720 directions, against the 280 mm the plate
    // has — 30.7 before #29 bounded the plate's own travel by the conditioning
    // of its Jacobian, which costs the loop a little authority at full tilt and
    // so lets the ball run half a millimetre further out.
    // This bound was raised to 90 mm when the ball first learned to hop,
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
// Where the tuning DOES start to let go is mapped elsewhere, because it is a
// different question: this sweep shoves a ball that has settled at the centre
// of a level plate, and the demo's Nudge button shoves one that is already
// tracking a circle with the legs already displaced.  That envelope, swept over
// speed as well as direction, is at `kMaxNudgeSpeed` — and at the bound it
// chooses the shipped tuning separates the ball in about one run in nine.
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
    // And with room to spare rather than by a whisker: 1.91 measured, a fifth
    // of a g.  The bar is 1.0, which still fails if the plate ever comes within
    // a tenth of a g of letting go at the disturbance the demo hands out.
    // Asserted only because `separated` is zero — see `Sweep::min_normal`.
    //
    // It was 1.62 before #29.  A plate that may not steer itself past a
    // Jacobian condition number of 20 also may not slam quite as hard, and
    // the margin the ball is held by is what that buys.
    ASSERT_TRUE(margin > 1.0);
}

// **What the servo rate limit does to this sweep, and — mostly — what it does
// not.**  #32 asks for the effect to be STATED under both tunings, so it is
// measured here rather than asserted in prose somewhere it can rot.
//
// The limit is a threshold on leg ERROR, not on ball speed: it engages above
// `rate_max * tau`, which against the shipped 10.47 rad/s servo and 0.05 s lag
// is 30.00 degrees.  Over 72 directions at the demo's own 0.26 m/s Nudge, the
// worst leg error any frame produced is
//
//     Nominal      27.30 deg   — 2.7 deg of headroom, the limit never engages
//     Aggressive   30.26 deg   — one frame past it, out of 187,920 leg-frames
//
// So the sweep is unchanged by it, on every number that has ever been asserted
// about it: nothing separates under either tuning, the normal-force margin is
// 1.913 m/s^2 Nominal and 1.354 Aggressive with and without the limit, and the
// peak reach is 31.15 mm and 28.57 mm either way.  The single Aggressive frame
// that ramps moves the worst leg rate in the sweep from 8.984 rad/s to 8.983.
//
// **Only the WITH-limit half of that is pinned below.**  The plant carries its
// rate limit as a property of the mechanism, so a test cannot step an unlimited
// plate without constructing a second one, and a `TableParams` built to be
// wrong is a fixture that outlives the question it answered.  The without-limit
// column was measured out of tree against `alpha_rate_max = infinity` and is
// reported here as history.  What IS pinned is the fact the whole comparison
// turns on: the worst leg error either tuning produces, against the handover.
//
// **That is the intended result, not a disappointment.**  Separation at this
// disturbance is an acceleration phenomenon and this is a rate limit; they bind
// in different regimes.  What the measurement buys is the knowledge that the
// limit sits *just* above where the Aggressive preset works, so it is insurance
// against the saturated recoveries this sweep does not contain rather than a
// tax on the demo the visitor actually sees.
void test_the_rate_limit_is_insurance_rather_than_a_tax() {
    const auto models = getBuiltinModels();
    const auto& e = cascadeModel(models);
    const TableParams tp = cascadeMechanism(e.params);
    // 30.00 degrees exactly, to within the float round-trip the model's tau
    // slider takes — the same micron-scale tolerance `samePlant` is sized by.
    const double handover = tp.alpha_rate_max * cascadeServoTau(e.params);
    ASSERT_NEAR(handover, 30.0 * M_PI / 180.0, 1e-6);

    struct Tuning { const char* name; Eigen::MatrixXd K; double worst_error; };
    Tuning tunings[] = {
        {"Nominal", defaultGain(e), 0.4765},
        {"Aggressive", gainForPreset(e, presetNamed("Aggressive")), 0.5282},
    };

    for (const Tuning& t : tunings) {
        double margin = 1e9;
        double leg_error = 0.0;
        int separated = 0;
        for (int i = 0; i < 72; ++i) {
            const Sweep sw = sweepOneDirection(e, t.K, Disturbance{},
                                               i * 2.0 * M_PI / 72.0);
            if (sw.separated) ++separated;
            margin = std::min(margin, sw.min_normal);
            leg_error = std::max(leg_error, sw.peak_leg_error);
        }
        ASSERT_EQ(separated, 0);
        ASSERT_TRUE(margin > 1.0);
        ASSERT_NEAR(leg_error, t.worst_error, 1e-3);
    }

    // Nominal stays clear of the limit; Aggressive touches it.  Pinned as the
    // ordering rather than only as two numbers, because it is the ordering
    // that says which tuning the limit is for.
    ASSERT_TRUE(tunings[0].worst_error < handover);
    ASSERT_TRUE(tunings[1].worst_error > handover);
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
//
// **RE-MEASURED AFTER #29, AND EVERY FAILURE IN THAT TABLE IS GONE.**  On the
// same 36 x 24 grid, with the plate no longer able to steer itself past a
// Jacobian condition number of 20:
//
//     speed   Nominal   Aggressive   Detuned      (failures out of 864)
//     0.25       0           0          0
//     0.30       0           0          0
//     0.35       0           0          0
//     0.40       0           0          0
//     0.45       0           0          0
//     0.50       0           0          0
//
// So every one of those losses was the plate changing assembly mode through a
// singularity rather than the gain running out of authority — the same defect,
// on the same grid, that #29's own sweep was failing on.  Both tables are kept
// because the OLD one is what chose `kMaxNudgeSpeed`, and a bound whose reason
// has moved is worth noticing rather than quietly adjusting.
//
// **AND RE-MEASURED AGAIN, AFTER THE LANDING BECAME A BOUNCE, WHICH IS WHERE
// THE BOUND ACTUALLY MOVED.**  #23's round two gives the landing a coefficient
// of restitution, and restitution moves the separation behaviour this whole
// envelope is made of: a shove hard enough to separate the ball used to end
// with it arriving and sticking, and now starts a train running about
// `e/(1-e) = 16` times the first flight.  For that second or so the ball is
// barely steerable — no rolling friction acts on it, and the plate can reach it
// only through the horizontal component of a normal impulse at each contact.
// Same 36 x 24 grid, against the contact model that ships:
//
//     speed   Nominal   Aggressive   Detuned      (failures out of 864)
//     0.18       0           0          0
//     0.20       0           0          0
//     0.22       0           1          0
//     0.24       0           0          0
//     0.25       4           2          0
//     0.30      10          35          0
//
// Detuned's row is again the reminder that these are slivers: 0.22 loses one
// direction and 0.24 loses none, which is what a 24-direction grid looks like
// walking past the side of one.  So the bound is the last speed BELOW the first
// failure of any tuning rather than the last clean row — **0.20**, by the same
// reading that put it at 0.30 against the table above.
//
// `kMaxSetpointSpeed` was re-measured in the same pass and did NOT move; see
// its own header for what it turned out to be about.  Both were the decision
// record's D16, and the ticket closing #23 owns them.
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

// **The no-pumping property, with the loop closed.**  This is what #23 closes
// on, and it is a claim about one scalar.
//
// Restitution acts on the RELATIVE normal velocity at the contact point, so a
// ball arriving at `u` onto a contact point rising at `p` leaves at
//
//     u_rebound = e u + p (1 + e)
//
// `p <= 0` therefore gives `u_rebound <= e u` — the bounce can never return
// more than it took — and a rising plate always adds energy, at 1.94x its own
// speed at `e = 0.94`.  That is a hard result about the impact law rather than
// a tuning, which is why the loop is handed a constraint (`holdContactDown`)
// instead of a fourth airborne target to chase.  The law itself is pinned in
// `test_ball_contact`, on a plate built by hand where the arithmetic is exact;
// what this pins is that the constraint SURVIVES the closed loop.
//
// **It is not asserted as an apex ratio, and the decision record's D10 asked
// for one.**  `e^2` bounds the ratio of successive apex HEIGHTS only over a
// plate that holds still, which is how the bounce was validated in
// `test_apexes_decay_at_e_squared_on_a_still_plate`.  With the loop running the
// plate does not hold still, and D8 is the reason: forbidden to rise, it
// descends under a ball it has just let go of, so the gap grows with no energy
// entering the ball at all.  Measured that way this sweep reports ratios near 2
// with `p` pinned to a thousandth of a metre per second — the metric is
// confounded by the very decision it sits beside.  So the mechanism is asserted
// directly, and D10's own secondary measures come with it as UPPER bounds:
// a loop that pumped would show it in bounce count and airborne frames, and
// those are pinned both ways rather than only used to prove the sweep is live.
//
// **Where the guarantee is given up, and by how much.**  `holdContactDown` asks
// for the leg rate that cancels the rise and stops there; two things can keep
// the plate from delivering it, and both are saturations rather than errors.
// The servo rate limit (#32) caps how fast a leg may move, and the retreat into
// the holdable set (#29) scales the whole triple back toward a level pose that
// sits higher than the one being asked for.  Measured over this sweep, frames
// where the contact point rose faster than 10 mm/s, out of every frame the ball
// entered airborne:
//
//     Nominal        0 of 2344     worst +0.0008 m/s
//     Aggressive    70 of 5359     worst +0.2959 m/s
//     Detuned        0 of 1871     worst +0.0009 m/s
//
// Aggressive is the tuning that slams the legs hardest, so it is the one that
// runs the servo out of speed, and the allowance below is written per tuning
// rather than as one loose number that would let Nominal rot quietly.
void test_the_loop_never_pumps_a_bouncing_ball() {
    const auto models = getBuiltinModels();
    const auto& e = cascadeModel(models);
    const SimPlate plate = cascadePlate(e.params);

    // **The allowances are pinned ON the measurement, not widened past it.**
    // Aggressive's 70 frames and 0.296 m/s are the servo's rate limit binding,
    // and a bound at twice either number would be an exemption rather than a
    // pin — it would let the residual double before anything failed.  10% of
    // headroom is for the arithmetic.
    struct Tuning {
        const char* name;
        Eigen::MatrixXd K;
        int allowed_rising;    ///< frames the contact point may rise on
        double worst_rise;     ///< [m/s] and by no more than this
        int max_bounces;       ///< and the train stays this short
        int max_airborne;
    };
    Tuning tunings[] = {
        {"Nominal", defaultGain(e), 0, 0.01, 700, 2600},
        {"Aggressive", gainForPreset(e, presetNamed("Aggressive")), 78, 0.33,
         1400, 5900},
        {"Detuned", gainForPreset(e, presetNamed("Detuned")), 0, 0.01, 620, 2100},
    };

    const double home_z = plate.kinematics()
                              .home_pose(cascadeHomeLegAngle(e.params)).z_c;

    for (Tuning& t : tunings) {
        int separated = 0, airborne = 0, rising = 0, impacts = 0, unended = 0;
        double worst_rise = 0.0, lowest_plate = home_z, end_of_run_z = home_z;
        for (int i = 0; i < 12; ++i) {
            for (int j = 0; j < 12; ++j) {
                SimInput in;
                in.design = cascadeDesign(e.params);
                in.design.K = t.K;
                in.path = openingPath();
                SimState s = simStart(plate, in.design.home_leg_rad,
                                      attractStart(in.path));
                const int at =
                    static_cast<int>((i / 12.0) * in.path.period_s / in.dt);
                const int total = at + static_cast<int>(20.0 / in.dt);
                const double theta = j * M_PI / 6.0;

                bool was_airborne = false, ever = false;
                for (int k = 0; k < total; ++k) {
                    if (k == at && !s.ball.airborne) {
                        s.ball.rolling(2) += kMaxNudgeSpeed * std::cos(theta);
                        s.ball.rolling(3) += kMaxNudgeSpeed * std::sin(theta);
                    }
                    const SimReport f = stepSim(plate, in, s);

                    // Only frames the ball ENTERED airborne are the
                    // constraint's.  The frame a separation happens on is one
                    // where the ball was still on the plate when the command
                    // was chosen — and the plate may well have been rising
                    // then, since it separated the ball by accelerating away
                    // from it rather than by descending.
                    if (was_airborne) {
                        if (f.contact_normal_rate > 0.01) ++rising;
                        worst_rise = std::max(worst_rise, f.contact_normal_rate);
                    }
                    if (f.impact_approach < 0.0) ++impacts;
                    if (f.airborne) {
                        ever = true;
                        ++airborne;
                        // **D8's recorded risk, measured rather than assumed.**
                        // "The plate descends at separation BECAUSE it
                        // accelerated away from the ball, so forbidding it to
                        // rise for a ~1 s train may strand it low and tilted."
                        lowest_plate = std::min(lowest_plate, s.pose.z_c);
                    }
                    was_airborne = f.airborne;
                    if (k + 1 == total) end_of_run_z = s.pose.z_c;

                    // The ball stays on the plate throughout, which is the
                    // claim `kMaxNudgeSpeed` is chosen to make.
                    ASSERT_TRUE(!f.left_plate);
                    if (k + 1 == total && f.airborne) ++unended;
                }
                if (ever) ++separated;
            }
        }
        // The sweep has to be measuring something.  A change that stopped the
        // ball separating at all would otherwise pass this test perfectly by
        // never testing it — which is how the 30-of-72 hop rate came to be
        // prose in the first place.
        ASSERT_TRUE(separated > 5);
        ASSERT_TRUE(impacts > 100);
        // And every train ends.  `bounceFloorSpeed` is what makes that true on
        // a still plate; this says the loop does not keep one alive.
        ASSERT_EQ(unended, 0);

        ASSERT_TRUE(rising <= t.allowed_rising);
        ASSERT_TRUE(worst_rise < t.worst_rise);

        // And the symptoms, as ceilings.  A loop that fed a train would show it
        // here first — more impacts, and longer in the air — so these fail
        // loudly rather than merely drifting.  Measured: 630 / 2344 under
        // Nominal, 1251 / 5359 under Aggressive, 548 / 1871 under Detuned.
        ASSERT_TRUE(airborne <= t.max_airborne);
        ASSERT_TRUE(impacts <= t.max_bounces);

        // **The risk D8 recorded did not bite.**  The plate drops 30.6 mm under
        // Nominal and 45.1 mm under Aggressive from a 212.1 mm home height, and
        // it is back at that height at the end of every run — so "stranded low"
        // is a dip of a fifth of the travel, not a floor the loop cannot climb
        // off.  Both halves are asserted: a plate that ran out of room and a
        // plate that never came back are different failures.
        ASSERT_TRUE(lowest_plate > home_z - 0.060);
        ASSERT_NEAR(end_of_run_z, home_z, 1e-3);
    }
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
// > **This was #29, and neither of the two answers it expected was right.**
// > It failed in direction 33 of 180 — 66.0 degrees — rolling the ball to
// > 293.7 mm against a rim at 280 mm.  The suspects were the analytic normal
// > force (#23) and `Q = 150` itself (#19), and it was neither: the plate came
// > out of a frame at Jacobian condition number 290 **on the other assembly**,
// > 13.8 degrees over with its downhill pointing at the ball, and spent the
// > next half second accelerating the ball it was trying to catch.  The loop
// > was asking for the right tilt the whole way and being handed its mirror.
// > `retreatToHoldable` is now bounded by `can_hold` rather than
// > `can_assemble`, the worst condition number this sweep reaches is 11.1, and
// > every direction holds.  See #29 and `onBuiltAssembly`.
void test_the_aggressive_preset_recovers_from_every_direction() {
    const auto models = getBuiltinModels();
    const auto& e = cascadeModel(models);

    const Eigen::MatrixXd K = gainForPreset(e, presetNamed("Aggressive"));

    for (int i = 0; i < 180; ++i)
        sweepOneDirection(e, K, Disturbance{}, i * M_PI / 90.0);
}

// And the fact that only exists now the ball can leave: an over-aggressive
// tuning does not merely overshoot, it takes the ball OFF THE SURFACE.  Q at
// 2000 slams the legs hard enough that the plate is moving out from under the
// ball rather than tilting under it.
//
// Measured against the application's own step, at the demo's own 0.26 m/s
// disturbance: the ball leaves the surface in 1 of 90 directions and rises
// 2.1 mm — and in 16 of 720, which is the sweep that says one is a real sliver
// rather than a grid artefact.  Before #23 the same tuning looked merely fast
// — 0.35 s to settle — because a ball glued to the plate cannot be thrown off
// it.
//
// **"And then off the plate" did not survive #29, and was never the gain's
// doing.**  This test used to assert a loss as well, and got one: the ball left
// the plate in 1 direction of 90.  That loss was the plate changing assembly
// mode through a singularity and tilting the ball away — the same defect, under
// a different gain, that the aggressive preset's sweep above was failing on.
// With the plate bounded to poses whose Jacobian it is allowed to believe, Q at
// 2000 keeps the ball at every direction of this disturbance.  That is asserted
// rather than dropped: it is the guard's own claim, and an absurd gain is where
// it is worth making.
//
// **The size of the hop is nothing like what this test used to claim.**  It
// said "one launch reached 689 mm of altitude".  Two of the reasons are in
// `test_the_shipped_tuning_holds_the_ball_through_its_own_disturbance` above —
// a differenced acceleration answering a stepping command with a delta
// function (#23), and a harness taking the servo rate at the wrong instant
// (#30).  A metre of altitude off a 0.26 m/s nudge was never physics.
//
// What survives, and is the point of the test, is that letting go of the ball
// at all is a thing an over-aggressive gain does and a sensible one does not:
// the three tunings the UI offers separate it in none of 720 directions at
// this same disturbance.  That is still the honest ceiling on #19's aggressive
// preset, and still a far better demonstration of over-aggressive control than
// a settling time.
//
// **Both of this test's own numbers moved when the landing became a bounce,
// and the interesting half is the direction count.**  It said the ball leaves
// in 1 of 90 directions and rises 2.1 mm; it now leaves in 90 of 90 and rises
// 28.7 mm.  The count did not move because the gain got worse — it moved
// because a separation that lasted one frame became one that lasts seconds,
// and a one-frame separation was invisible to a test that samples
// `SimReport::airborne` at the frame boundary.  The ball was leaving all along;
// under `kRestitution` it stays gone long enough to be seen.
//
// Which also relocates the cause.  This sweep settles the ball from 60, -40 mm
// before it shoves it, and at Q = 2000 the ball separates on frame 11 of that
// settling — the same frame in every direction, because the shove has not
// happened yet.  So what the 90 counts is the tuning throwing the ball while
// merely *centring* it, which is a stronger statement about the gain than the
// disturbance sweep was making.
//
// The ball still comes back and stays on the plate in every direction, and the
// plate still never leaves the assembly it is built in.
void test_an_over_aggressive_tuning_takes_the_ball_off_the_surface() {
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

    int lost = 0, separated = 0;
    double peak_hop = 0.0;
    for (int i = 0; i < 90; ++i) {
        const double theta = i * M_PI / 45.0;
        SimState st = simStart(plate, in.design.home_leg_rad,
                               Eigen::Vector4d(s.start_x, s.start_y, 0.0, 0.0));

        bool airborne_here = false;
        for (int k = 0; k < settle + static_cast<int>(6.0 / dt); ++k) {
            if (k == settle && !st.ball.airborne) {
                st.ball.rolling(2) += s.speed * std::cos(theta);
                st.ball.rolling(3) += s.speed * std::sin(theta);
            }
            const SimReport frame = stepSim(plate, in, st);
            // Four times the aggressive preset is where the plate is most
            // likely to be steered somewhere it should not be, so it is the
            // best place to ask: it is still the assembly the machine is built
            // in.  See #29 and `onBuiltAssembly`.
            ASSERT_TRUE(onBuiltAssembly(plate.kinematics(), st.alpha_rad, st.pose));
            if (frame.airborne) airborne_here = true;
            peak_hop = std::max(peak_hop, frame.ball_plate(2) - kBallRadius);
            if (frame.left_plate) { ++lost; break; }
        }
        if (airborne_here) ++separated;
    }
    // It really does take the ball off the surface — every one of these 90
    // directions, because the separation is in the settling they share rather
    // than in the shove that distinguishes them.  Asserted as "all of them"
    // rather than "at least one": a tuning that only threw the ball on a
    // shove would be a materially better tuning than this one, and should
    // fail here rather than pass quietly.
    ASSERT_EQ(separated, 90);
    // And does not lose it, at the disturbance the demo's own buttons offer.
    ASSERT_EQ(lost, 0);
    // Millimetres, not metres: 28.7 mm measured over these 90 directions,
    // against the 689 mm this test once claimed and the 2.1 mm it claimed
    // while the ball could not bounce.  Both bounds are assertions — too small
    // and the plate has stopped letting go at all, too large and something is
    // differencing a step again.
    ASSERT_TRUE(peak_hop > 0.005);
    ASSERT_TRUE(peak_hop < 0.060);
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
    test_the_rate_limit_is_insurance_rather_than_a_tax();
    test_a_shove_while_tracking_is_rejected_from_every_direction();
    test_the_nudge_buttons_cannot_compose_past_the_bound();
    test_the_loop_never_pumps_a_bouncing_ball();
    test_the_aggressive_preset_recovers_from_every_direction();
    test_an_over_aggressive_tuning_takes_the_ball_off_the_surface();
    test_ten_minutes_unattended_never_loses_the_ball();
    std::printf("test_attract_mode: all passed\n");
    return 0;
}
