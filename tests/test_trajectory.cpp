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
#include "sim_step.h"
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

/// `stepSim` — the step the application drives, with the setpoint driven by a
/// path and the loop given the path's velocity as its reference velocity.
///
/// `feedforward` chooses whether that velocity is SUPPLIED, and no longer how
/// it is applied: the arithmetic lives in `legCommand`, where the application
/// gets it from too.  Withholding it is the one thing this harness does that
/// the application never does, so it is expressed the only way it can be — by
/// telling the step there is no path and holding the setpoint at the point the
/// path is passing through this frame.  A held setpoint has no velocity, which
/// is exactly the case being measured.
///
/// Nothing else here is a copy of anything.  This file used to carry the whole
/// causal order, and a second implementation of the loop is what let a real bug
/// hide once already (#22) — a harness that recomputes the thing it is meant to
/// be measuring is not evidence about the application.  See #30.
Track follow(const ModelEntry& e, const Eigen::MatrixXd& K,
             const SetpointPath& asked, bool feedforward, double duration) {
    const SimPlate plate = cascadePlate(e.params);

    // The fillet the step is going to stamp on anyway, taken here as well so
    // that the two things this harness reads the path for itself — where to
    // start the ball, and where the held setpoint is each frame — are on the
    // same curve the loop is being driven round.  Reading it off the plate
    // rather than copying a number is the point: a harness with its own idea
    // of `accel_max` would be measuring a path the application never runs.
    const SetpointPath path = plate.feasible(asked);

    SimInput in;
    in.design = cascadeDesign(e.params);
    in.design.K = K;
    in.path = path;

    const Eigen::Vector2d start = pathPoint(path, 0.0);
    SimState s = simStart(plate, in.design.home_leg_rad,
                          Eigen::Vector4d(start(0), start(1), 0.0, 0.0));

    const double dt = in.dt;
    Track r;
    double t = 0.0, sum = 0.0;
    int n = 0;
    for (int k = 0; k < static_cast<int>(duration / dt); ++k) {
        if (!feedforward) {
            // A HELD setpoint, moved by hand to where the path is passing this
            // frame.  Held is the whole trick: a held setpoint has no velocity,
            // which is exactly the case being measured, and the loop is handed
            // the same position it would have been either way.
            //
            // The phase has to be walked here, because a `Fixed` path is one
            // `stepSim` does not advance — nothing is running a lap.  That is
            // the harness reaching into the step's state, and it is worth being
            // plain that it is: what stops it becoming a second copy of the
            // timing rule is that it walks the phase with `advancePhase`, the
            // function `stepPath` itself walks it with, at the same `dt`.  The
            // rule has one implementation; this drives it directly instead of
            // through a path that is not running.
            in.path.shape = PathShape::Fixed;
            in.held_setpoint = pathPoint(path, s.path_phase);
            s.path_phase = advancePhase(s.path_phase, dt, path.period_s);
        }
        const SimReport frame = stepSim(plate, in, s);
        const Eigen::Vector2d sp =
            feedforward ? frame.setpoint : in.held_setpoint;
        t += dt;

        r.max_radius = std::max(r.max_radius,
                                std::hypot(frame.ball_plate(0), frame.ball_plate(1)));
        if (t > path.period_s) {          // after one lap, so the start is not counted
            const double err = std::hypot(frame.ball_plate(0) - sp(0),
                                          frame.ball_plate(1) - sp(1));
            r.max_error = std::max(r.max_error, err);
            sum += err;
            ++n;
        }
        if (frame.left_plate) {
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
    // The lap floor moves with the fillet — a path that rounds its corners has
    // less ground to cover — so the floor has to be asked of a path that knows
    // what the plate can do.  See `minPeriod`.
    const double a_max = cascadePlate(e.params).maxBallAccel();

    double worst_radius = 0.0, worst_error = 0.0;
    for (PathShape s : {PathShape::Circle, PathShape::Square, PathShape::Triangle}) {
        for (double r_mm : {20.0, 60.0, 120.0, kMaxPathRadius * 1000.0}) {
            SetpointPath p;
            p.shape = s;
            p.radius_m = r_mm * 1e-3;
            p.accel_max = a_max;
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
    // Measured, all 24 settings kept: the ball reaches 192 mm on the largest
    // path, against the 280 mm the plate has, and its worst mean error is
    // 21 mm — at the fastest triangle, which is the setting that is SUPPOSED
    // to look hard.
    //
    // Both numbers came down when the corners were filleted (#31).  The mean
    // error was 36 mm and is 21; the reach was 195 mm and is 192.  Neither is
    // the point of the fillet and both are consequences of it: the reference
    // stopped asking for a turn the ball cannot make, so the ball stopped
    // being thrown wide of one.
    ASSERT_TRUE(worst_radius < 0.24);
    ASSERT_TRUE(worst_error < 0.03);
}

// What the corner actually hands the actuator, which is the whole of #31.
//
// A sharp corner turns the reference velocity through a right angle between
// two frames.  Sampled at 60 Hz on the fastest square that is 0.354 m/s in one
// frame — an implied 21 m/s^2, against the 1.89 the plate can give the ball —
// and the reference is asking for infinity, the sampling is merely what stops
// the number being one.  That impulse reaches the legs through K's velocity
// columns whether or not the ball can follow it.
//
// Blended, the reference's own acceleration is `v^2/rho = accel_max` through
// the corner and zero along the straights, by construction.  This measures it
// where it matters — off `SimReport::setpoint_velocity`, the quantity the loop
// was actually handed, rather than off the path re-read afterwards.
void test_the_reference_never_asks_for_more_than_the_ball_can_do() {
    const auto models = getBuiltinModels();
    const auto& e = cascadeModel(models);
    const SimPlate plate = cascadePlate(e.params);
    const double a_max = plate.maxBallAccel();

    SimInput in;
    in.design = cascadeDesign(e.params);
    in.design.K = defaultGain(e);

    for (PathShape shape : {PathShape::Square, PathShape::Triangle}) {
        SetpointPath p;
        p.shape = shape;
        p.radius_m = kMaxPathRadius;
        p.accel_max = a_max;
        p.period_s = minPeriod(p);          // the fastest the sliders offer
        in.path = p;

        SimState st = simStart(plate, in.design.home_leg_rad);
        Eigen::Vector2d previous = Eigen::Vector2d::Zero();
        double worst = 0.0;
        const int frames = static_cast<int>(2.0 * p.period_s / in.dt);
        for (int k = 0; k < frames; ++k) {
            const SimReport f = stepSim(plate, in, st);
            if (k > 0)
                worst = std::max(worst,
                                 (f.setpoint_velocity - previous).norm() / in.dt);
            previous = f.setpoint_velocity;
        }
        // A frame straddling the join between a fillet and a straight averages
        // the two, so a 60 Hz difference can only ever UNDER-read the peak —
        // it cannot manufacture one.  Measured 1.87 m/s^2 against the 1.888
        // the plate can give.  The slack is for the arithmetic, not for a
        // margin: 20% here would still be six times under the sharp corner.
        ASSERT_TRUE(worst < 1.2 * a_max);
    }
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
// than a circle of the same size and lap time, and it has to be — a corner
// turns much harder than a circle of the same period does, and a
// bandwidth-limited loop rounds what it cannot turn into.
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
    // Measured: the square's worst error is 7.8 mm against the circle's 5.1 mm
    // on the same size and the same lap time.  Half again, and it is the
    // corner that does it.
    //
    // **This was 14.2 mm and nearly three times, before the corners were
    // filleted.**  The 3.7 mm fillet a 120 mm square at an eight-second lap
    // takes is small enough that the shape still reads as a square and the
    // corner is still visibly the hard part — the demo keeps its point — but
    // half of what the ball used to lose there was the reference asking for a
    // turn no tilt of this plate could produce.  See #31.
    //
    // The assertion is a RATIO rather than either number, because the square's
    // figure is a peak sampled at 60 Hz near the tightest part of the path.
    // The fillet makes it far less frame-sensitive than the step did — the old
    // number swung between 14.2 and 14.9 mm on a phase perturbation of one bit
    // — but it is still a peak, and a fourth significant figure on it would be
    // a claim about frame alignment rather than about the corner.
    ASSERT_TRUE(q.max_error > 1.3 * c.max_error);
    // And the fillet really did take the bulk of it: a sharp corner put the
    // square nearly three times the circle, so anything at or above that is
    // the fillet having been lost.
    ASSERT_TRUE(q.max_error < 2.2 * c.max_error);
}

}  // namespace

int main() {
    test_every_offered_setting_keeps_the_ball();
    test_the_reference_never_asks_for_more_than_the_ball_can_do();
    test_feedforward_is_worth_having();
    test_a_corner_is_harder_than_a_curve();
    std::printf("test_trajectory: all passed\n");
    return 0;
}
