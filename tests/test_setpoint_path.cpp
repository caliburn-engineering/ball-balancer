// tests/test_setpoint_path.cpp
//
// The paths the ball is asked to trace.  Nothing here closes a loop — these are
// claims about the shapes themselves, and they are worth making separately
// because a tracking failure is much easier to read when the thing being
// tracked is known to be right.
//
// The last section is the exception, and deliberately so: `accel_max` is the
// one thing about a path that comes from the PLANT, and a file that only ever
// made up its own value could not catch the fillet being sized off a plate
// nobody ships.  See #31.
#include "setpoint_path.h"

#include "cascade_fixture.h"
#include "sim_step.h"
#include "table_kinematics.h"
#include "trajectory_controls.h"
#include "test_helpers.h"

#include <algorithm>
#include <cmath>
#include <limits>

using namespace caliburn;

namespace {

/// A round number standing in for what the plate can give the ball — the
/// shipped geometry answers 1.888 m/s^2, and the arithmetic below is about the
/// fillet rather than about that plate.  The last section pins the real one.
constexpr double kSomeAccel = 2.0;

SetpointPath shape(PathShape s, double r = 0.12, double period = 10.0,
                   double accel_max = 0.0) {
    SetpointPath p;
    p.shape = s;
    p.radius_m = r;
    p.period_s = period;
    p.accel_max = accel_max;
    return p;
}

// Every shape closes: a lap returns exactly where it started, at every size
// and period.  A path that drifts would send the ball off the plate eventually
// and the drift would be invisible for the first few laps.
void test_every_path_closes() {
    for (PathShape s : {PathShape::Circle, PathShape::Square, PathShape::Triangle}) {
        for (double r : {0.04, 0.12, 0.20}) {
            for (double T : {3.0, 10.0, 30.0}) {
                const SetpointPath p = shape(s, r, T);
                for (int lap = 1; lap <= 4; ++lap) {
                    const Eigen::Vector2d a = pathPoint(p, 0.7);
                    const Eigen::Vector2d b = pathPoint(p, 0.7 + lap);
                    ASSERT_NEAR((a - b).norm(), 0.0, 1e-12);
                }
            }
        }
    }
}

// Nothing ever goes outside the radius asked for.  This is what lets the plate
// bound be set once, against the radius slider, rather than per shape.
void test_no_path_leaves_its_radius() {
    for (PathShape s : {PathShape::Circle, PathShape::Square, PathShape::Triangle}) {
        const SetpointPath p = shape(s, 0.12);
        for (int i = 0; i < 2000; ++i)
            ASSERT_TRUE(pathPoint(p, i / 2000.0).norm() <= 0.12 + 1e-12);
    }
}

// And each one actually reaches it, at its corners.  Sizing from the
// circumradius is what makes the shapes comparable at the same setting.
void test_every_shape_reaches_its_radius() {
    for (PathShape s : {PathShape::Circle, PathShape::Square, PathShape::Triangle}) {
        const SetpointPath p = shape(s, 0.12);
        double furthest = 0.0;
        for (int i = 0; i < 2000; ++i)
            furthest = std::max(furthest, pathPoint(p, i / 2000.0).norm());
        ASSERT_NEAR(furthest, 0.12, 1e-6);
    }
}

// A filleted polygon reaches slightly less far, because the fillet cuts the
// corner off — and by exactly `rho*(sec(pi/n) - 1)`, which is the corner's
// distance from its own fillet circle.  Worth pinning rather than bounding:
// the size slider is still measured from the corner, and how much of the
// corner the fillet eats is the number that says whether that is honest.
void test_a_filleted_polygon_falls_short_of_its_radius_by_the_fillet() {
    for (PathShape s : {PathShape::Square, PathShape::Triangle}) {
        const int n = pathCorners(s);
        for (double T : {0.0, 10.0}) {
            SetpointPath p = shape(s, kMaxPathRadius, 10.0, kSomeAccel);
            p.period_s = (T > 0.0) ? T : minPeriod(p);   // the floor, and slow
            const double expect =
                p.radius_m - filletRadius(p) * (1.0 / std::cos(M_PI / n) - 1.0);
            double furthest = 0.0;
            for (int i = 0; i < 20000; ++i)
                furthest = std::max(furthest, pathPoint(p, i / 20000.0).norm());
            ASSERT_NEAR(furthest, expect, 1e-6);
            // Measured on the 180 mm shapes: the square gives up 12.9 mm at
            // its floor and 2.1 mm at a ten-second lap, and the triangle
            // 31.3 mm and 4.2 mm — a triangle's corner is sharper, so its
            // fillet cuts deeper into it.  `sec(pi/3) - 1` is exactly 1, which
            // is why a triangle's shortfall IS its fillet radius.
            ASSERT_TRUE(p.radius_m - furthest < 0.032);
        }
    }
}

// Constant speed, not constant angle.  A polygon traversed by sweeping an
// angle would crawl through its corners, which is precisely where the
// interesting behaviour is.
void test_the_polygons_are_walked_at_constant_speed() {
    for (PathShape s : {PathShape::Square, PathShape::Triangle}) {
        // Blended and sharp alike: the fillet changes the direction the
        // setpoint turns through, not the rate it covers ground at.  The
        // filleted sweep can include the corners, because the filleted path has
        // a velocity there.
        for (double a : {0.0, kSomeAccel}) {
            const SetpointPath p = shape(s, 0.12, 6.0, a);
            const double expect = pathLength(p) / p.period_s;
            for (int i = 0; i < 600; ++i) {
                const double u = (a > 0.0) ? i / 600.0 : (i + 0.5) / 600.0;
                ASSERT_NEAR(pathVelocity(p, u).norm(), expect, 1e-9);
            }
        }
    }
}

// The velocity is the position's derivative with respect to TIME, which is the
// claim the feedforward rests on.  The position is a function of phase and the
// phase advances at `1 / period_s`, so the chain rule is the thing being
// checked: `dp/dt = (dp/du) / period_s`.  Away from the corners, where the
// derivative genuinely does not exist.
void test_velocity_is_the_derivative_of_position() {
    const double h = 1e-6;
    for (PathShape s : {PathShape::Circle, PathShape::Square, PathShape::Triangle}) {
        const SetpointPath p = shape(s, 0.12, 7.0);
        for (int i = 0; i < 200; ++i) {
            const double u = (i + 0.5) / 200.0;
            const Eigen::Vector2d fd =
                (pathPoint(p, u + h) - pathPoint(p, u - h)) / (2 * h * p.period_s);
            ASSERT_TRUE((fd - pathVelocity(p, u)).norm() < 1e-4);
        }
    }
}

// And on a filleted path it holds AT the corners too, which is the whole reason
// the geometry is filleted rather than the velocity slewed.  Slewing
// `pathVelocity` alone would have made this test fail at every corner — the
// reference velocity would have stopped being the reference position's
// derivative, which is a new lie in exactly the place #31 is removing one.
void test_the_fillet_keeps_the_velocity_the_position_s_derivative_at_the_corners() {
    const double h = 1e-7;
    for (PathShape s : {PathShape::Square, PathShape::Triangle}) {
        const int n = pathCorners(s);
        const SetpointPath p = shape(s, 0.12, 4.0, kSomeAccel);
        // Straight through each corner's phase, and either side of the joins
        // where a fillet meets a straight.
        for (int k = 0; k < n; ++k) {
            for (double off : {-1e-3, -1e-5, 0.0, 1e-5, 1e-3}) {
                const double u = k / static_cast<double>(n) + off;
                const Eigen::Vector2d fd =
                    (pathPoint(p, u + h) - pathPoint(p, u - h)) / (2 * h * p.period_s);
                ASSERT_TRUE((fd - pathVelocity(p, u)).norm() < 1e-4);
            }
        }
    }
}

// A SHARP corner is a step in the reference velocity — that is what a corner
// means.  Kept, because it is the thing the fillet is measured against: #31
// reversed the claim that the step was harmless, not the claim that it is
// there.  `accel_max = 0` is the reference this code drove before the fillet.
void test_a_sharp_corner_is_a_step_in_the_reference_velocity() {
    const SetpointPath p = shape(PathShape::Square, 0.12, 8.0);
    const double corner_u = 0.25;                  // end of the first edge
    const Eigen::Vector2d before = pathVelocity(p, corner_u - 1e-5);
    const Eigen::Vector2d after = pathVelocity(p, corner_u + 1e-5);
    // A square's corner turns the velocity through a right angle.
    ASSERT_NEAR(before.normalized().dot(after.normalized()), 0.0, 1e-6);
    ASSERT_NEAR(before.norm(), after.norm(), 1e-9);
}

// And a filleted one is not.  The reference velocity turns through the same
// right angle over the fillet instead of between two samples of it, so the
// direction moves continuously and the speed never changes at all.
void test_a_filleted_corner_turns_the_reference_velocity_continuously() {
    for (PathShape s : {PathShape::Square, PathShape::Triangle}) {
        const int n = pathCorners(s);
        const SetpointPath p = shape(s, 0.12, 4.0, kSomeAccel);
        const double speed = pathLength(p) / p.period_s;
        for (int k = 0; k < n; ++k) {
            const double u = k / static_cast<double>(n);
            const Eigen::Vector2d before = pathVelocity(p, u - 1e-6);
            const Eigen::Vector2d after = pathVelocity(p, u + 1e-6);
            // Two microphases apart the reference velocity has barely moved —
            // measured, a ten-thousandth of the speed it is travelling at.
            // The sharp corner below moves it by 1.41 times that speed across
            // the same gap, more than a thousand times as much, and that
            // difference is the whole of the ticket.
            ASSERT_TRUE((after - before).norm() < 1e-3 * speed);
            // And the speed itself does not change at all.
            ASSERT_NEAR(before.norm(), after.norm(), 1e-12);
        }
        SetpointPath sharp = p;
        sharp.accel_max = 0.0;
        ASSERT_TRUE((pathVelocity(sharp, 1e-6) - pathVelocity(sharp, -1e-6)).norm() >
                    1.0 * pathLength(sharp) / sharp.period_s);
        // The turn still happens — it is spread over the fillet, not removed.
        // Over one corner's fillet the velocity comes round by the full
        // exterior angle, 2*pi/n.
        const double half = 0.5 / n * (filletRadius(p) * 2.0 * M_PI / n) /
                            (pathLength(p) / n);
        const Eigen::Vector2d in_ = pathVelocity(p, -half);
        const Eigen::Vector2d out = pathVelocity(p, half);
        ASSERT_NEAR(in_.normalized().dot(out.normalized()),
                    std::cos(2.0 * M_PI / n), 1e-6);
    }
}

// The number the whole ticket is: `rho = v^2 / a_max`, at the speed the path is
// actually walked — which is the speed of the BLENDED path, since filleting it
// shortened it.  Solved rather than iterated, so this checks the fixed point
// holds rather than that an iteration converged.
void test_the_fillet_radius_is_v_squared_over_a_max() {
    for (PathShape s : {PathShape::Square, PathShape::Triangle}) {
        for (double r : {0.02, 0.06, 0.12, kMaxPathRadius}) {
            for (double T : {2.0, 5.0, 30.0}) {
                const SetpointPath p = shape(s, r, T, kSomeAccel);
                const double v = pathLength(p) / p.period_s;
                ASSERT_NEAR(filletRadius(p), v * v / kSomeAccel, 1e-12);
            }
        }
    }
}

// Which means it shrinks away on a slow lap and opens up on a fast one — #24's
// bandwidth argument made visible in the TARGET rather than inferred from the
// ball's overshoot.
void test_the_fillet_grows_as_the_lap_tightens() {
    const SetpointPath fast = shape(PathShape::Square, kMaxPathRadius, 4.0, kSomeAccel);
    const SetpointPath slow = shape(PathShape::Square, kMaxPathRadius, 30.0, kSomeAccel);
    // Measured: 29.3 mm at a four-second lap and 0.57 mm at thirty.  Under a
    // millimetre is a corner nobody can see, which is the right answer for a
    // lap slow enough that the ball has no trouble with one.
    ASSERT_TRUE(filletRadius(fast) > 0.028 && filletRadius(fast) < 0.033);
    ASSERT_TRUE(filletRadius(slow) < 0.001);
    // Roughly the square of the lap ratio, since the speed is roughly the
    // length over the period and the length barely moves.
    ASSERT_TRUE(filletRadius(fast) / filletRadius(slow) > 40.0);
}

// The reference's own acceleration, which is what "feasible" means here: the
// fillet turns at exactly `a_max` and the straights do not turn at all, so the
// most the reference ever asks of the ball is the most the plate can give it.
void test_the_filleted_reference_never_exceeds_a_max() {
    const double h = 1e-6;
    for (PathShape s : {PathShape::Square, PathShape::Triangle}) {
        const SetpointPath p = shape(s, kMaxPathRadius, 4.0, kSomeAccel);
        double worst = 0.0;
        for (int i = 0; i < 20000; ++i) {
            const double u = i / 20000.0;
            const Eigen::Vector2d a =
                (pathVelocity(p, u + h) - pathVelocity(p, u - h)) / (2 * h * p.period_s);
            worst = std::max(worst, a.norm());
        }
        // Exactly `a_max` on the fillets, to the difference's own resolution.
        ASSERT_NEAR(worst, kSomeAccel, 1e-4);
    }

    // And what it replaced.  A sharp corner turns the velocity through a right
    // angle in no distance at all, so the same difference reports whatever the
    // step size lets it — 1e5 m/s^2 at this one, and unbounded as `h` shrinks.
    // That is the modelling error, not a large number.
    const SetpointPath sharp = shape(PathShape::Square, kMaxPathRadius, 4.0);
    const Eigen::Vector2d step =
        (pathVelocity(sharp, 0.25 + h) - pathVelocity(sharp, 0.25 - h)) /
        (2 * h * sharp.period_s);
    ASSERT_TRUE(step.norm() > 1e4 * kSomeAccel);
}

// The fillet has a ceiling: the polygon's inradius, where the fillets meet each
// other and the shape has become its own incircle.  Nothing the sliders offer
// comes near it — the 180 mm square asks 33 mm against a 127 mm cap — but the
// cap is what makes an absurd `accel_max` degrade into a circle rather than
// into a shape that folds through itself.
void test_the_fillet_is_capped_at_the_incircle() {
    for (PathShape s : {PathShape::Square, PathShape::Triangle}) {
        const int n = pathCorners(s);
        const SetpointPath p = shape(s, 0.12, 10.0, 1e-6);   // an absurdly weak plate
        const double inradius = 0.12 * std::cos(M_PI / n);
        ASSERT_NEAR(filletRadius(p), inradius, 1e-12);
        // Fully filleted IS the incircle: constant radius, and a perimeter to
        // match.
        for (int i = 0; i < 1000; ++i)
            ASSERT_NEAR(pathPoint(p, i / 1000.0).norm(), inradius, 1e-9);
        ASSERT_NEAR(pathLength(p), 2.0 * M_PI * inradius, 1e-9);
    }
}

// Blending shortens the path, because the two tangent lengths a corner gives
// up are longer than the arc it gets back.  That is why the lap floor comes
// DOWN — a path with rounded corners has less ground to cover at the same cap.
void test_a_filleted_polygon_is_shorter_and_so_is_its_fastest_lap() {
    for (PathShape s : {PathShape::Square, PathShape::Triangle}) {
        const int n = pathCorners(s);
        SetpointPath sharp = shape(s, kMaxPathRadius, 4.0);
        SetpointPath filleted = sharp;
        filleted.accel_max = kSomeAccel;

        const double rho = filletRadius(filleted);
        const double saved = rho * (2.0 * n * std::tan(M_PI / n) - 2.0 * M_PI);
        ASSERT_NEAR(pathLength(sharp) - pathLength(filleted), saved, 1e-12);
        ASSERT_TRUE(pathLength(filleted) < pathLength(sharp));

        // And the floor: measured, the 180 mm square's fastest offered lap
        // moves from 4.07 s to 3.86 s.
        ASSERT_TRUE(minPeriod(filleted) < minPeriod(sharp));
        // The floor is a fixed point of a circular definition — the fillet sets
        // the length, the length sets the lap, the lap sets the fillet — so the
        // thing to check is that it closes: a path AT its floor runs at exactly
        // the cap.
        filleted.period_s = minPeriod(filleted);
        ASSERT_NEAR(pathLength(filleted) / filleted.period_s, kMaxSetpointSpeed, 1e-12);
    }
    // Except where `kMinLapSeconds` is the binding bound instead, which is the
    // small paths — there the floor is two seconds whatever the fillet does.
    SetpointPath small = shape(PathShape::Square, 0.02, 10.0, kSomeAccel);
    ASSERT_NEAR(minPeriod(small), kMinLapSeconds, 1e-12);
}

// A degenerate period must not divide, and `Fixed` must stay put — it is the
// setpoint the balance loop has always had, and the path layer has to leave it
// alone rather than dragging it to the origin every frame.
void test_the_degenerate_cases_behave() {
    SetpointPath p = shape(PathShape::Circle, 0.12, 0.0);
    ASSERT_TRUE(std::isfinite(pathPoint(p, 0.3).norm()));
    ASSERT_NEAR(pathVelocity(p, 0.3).norm(), 0.0, 1e-15);
    // And a zero lap does not advance the phase either, rather than dividing.
    ASSERT_NEAR(advancePhase(0.3, 1.0 / 60.0, 0.0), 0.3, 1e-15);

    p = shape(PathShape::Fixed);
    ASSERT_NEAR(pathPoint(p, 0.3).norm(), 0.0, 1e-15);
    ASSERT_NEAR(pathLength(p), 0.0, 1e-15);

    // A non-finite phase cannot come out of `advancePhase`, but `pathPoint` is
    // public and the polygon walk casts to `int` — undefined for a NaN rather
    // than merely wrong.  It is stopped at the door, so this is a real branch
    // and gets a real test.
    for (PathShape s : {PathShape::Circle, PathShape::Square, PathShape::Triangle})
        ASSERT_TRUE(std::isfinite(pathPoint(shape(s), NAN).norm()));
}

// ---------------------------------------------------------------------------
// The sliders, mid-run (#24 round two)
// ---------------------------------------------------------------------------
//
// The defect these exist for: the phase used to be `t / period_s`, derived
// from the whole time the simulation had been running.  Changing the lap moved
// it by `t dT / T^2`, so after 100 s a nudge from 10.0 to 9.5 s threw the
// setpoint 170 degrees round the path and the loop hauled the ball clean
// across the plate after it.  It grew with run time and it wrapped, so any
// nudge could land anywhere.
//
// Every test in this file used to evaluate the path at fixed parameters, which
// is exactly why not one of them saw it.  These change a slider mid-run.
//
// **They drive `TrajectoryControls`, which is the object the panel drives.**
// They used to drive a `Run` harness that re-stated the panel's three rules
// beside it — the lap floor, the re-seeded phase, the size slider moving the
// floor with it — which is a test of the harness's copy rather than of the
// application (#28).  The rules are the object's methods now, and what is left
// here is a clock.

/// Where a frame would have left the setpoint: on the path, at `phase`.
///
/// Through `fromReport`, because that is how the setpoint reaches the controls
/// in the application — the step reports where it drove the setpoint to and the
/// panel holds the answer.  The tests need it because `setShape` re-seeds from
/// where the setpoint IS, and that is the whole of what makes a shape change
/// not jump.
void frameLeftItAt(TrajectoryControls& tc, double phase) {
    SimReport r;
    r.setpoint = pathPoint(tc.path(), phase);
    r.setpoint_velocity = pathVelocity(tc.path(), phase);
    tc.fromReport(r);
}

/// Run the lap for `seconds` at the application's own frame rate.
///
/// The phase is walked with `advancePhase`, which is the function `stepPath`
/// itself walks it with at the same `dt`, rather than with a second expression
/// for the same thing.  Nothing else is here: every rule these tests are about
/// is a method on `TrajectoryControls` or a function in `setpoint_path.h`, and
/// what is left is a clock.
///
/// The plate is deliberately absent, because these tests choose `accel_max` —
/// a sharp-cornered path is where "the setpoint does not move" is EXACT, and a
/// filleted one changes shape when the lap changes.  The section after this one
/// drives `stepSim` on the shipped plate instead, which is where the wiring
/// claims belong.
void advance(TrajectoryControls& tc, double& phase, double seconds) {
    const double dt = 1.0 / 60.0;
    for (int i = 0; i < static_cast<int>(seconds / dt); ++i)
        phase = advancePhase(phase, dt, tc.path().period_s);
    frameLeftItAt(tc, phase);
}

// The lap slider changes the RATE and nothing else.  Not approximately —
// exactly: the phase is a number the slider does not touch.
void test_changing_the_lap_leaves_the_setpoint_where_it_is() {
    for (PathShape s : {PathShape::Circle, PathShape::Square, PathShape::Triangle}) {
        // Deliberately at four run lengths spanning four orders of magnitude.
        // The old failure was proportional to `t`, so a test that only ever
        // ran for a second would have passed against the bug.
        for (double run_s : {0.5, 10.0, 100.0, 1000.0}) {
            TrajectoryControls tc(shape(s, 0.12, 10.0));
            double phase = 0.0;
            advance(tc, phase, run_s);

            const Eigen::Vector2d before = pathPoint(tc.path(), phase);
            tc.setLapS(9.5);
            const Eigen::Vector2d after = pathPoint(tc.path(), phase);
            ASSERT_NEAR((after - before).norm(), 0.0, 1e-15);

            // And again the other way, and by a lot rather than a nudge.
            tc.setLapS(30.0);
            ASSERT_NEAR((pathPoint(tc.path(), phase) - after).norm(), 0.0, 1e-15);
        }
    }
}

// The reference velocity, on the other hand, IS allowed to step — and must.
// The setpoint really was just asked to travel at a different speed, and a
// feedforward that ignored that would be feeding forward the old lap.
void test_the_reference_velocity_steps_on_a_lap_change() {
    // Every shape: the polygons' velocity is just as period-dependent as the
    // circle's, and the position tests beside this one sweep all three.
    for (PathShape s : {PathShape::Circle, PathShape::Square, PathShape::Triangle}) {
        TrajectoryControls tc(shape(s, 0.12, 10.0));
        double phase = 0.0;
        advance(tc, phase, 100.0);

        const Eigen::Vector2d before = pathVelocity(tc.path(), phase);
        tc.setLapS(5.0);                    // above every shape's floor at 120 mm
        const Eigen::Vector2d after = pathVelocity(tc.path(), phase);
        ASSERT_NEAR(tc.lapS(), 5.0, 1e-12);

        // Same direction — the setpoint has not turned, it has sped up.
        ASSERT_NEAR(before.normalized().dot(after.normalized()), 1.0, 1e-12);
        // Twice the lap rate, twice the speed.
        ASSERT_NEAR(after.norm(), 2.0 * before.norm(), 1e-12);
    }
}

// The size slider slides the setpoint radially and does not rotate it: same
// angle, longer radius, which is the same point on a bigger shape.  Checked
// both ways round the floor — growing a path pushes the lap up with it, and
// that indirect change of `period_s` was the second route into the bug.
void test_changing_the_size_moves_the_setpoint_radially_only() {
    for (PathShape s : {PathShape::Circle, PathShape::Square, PathShape::Triangle}) {
        for (double run_s : {0.5, 100.0, 1000.0}) {
            TrajectoryControls tc(shape(s, 0.02, 10.0));
            tc.setLapS(minPeriod(tc.path()));     // sitting on the floor
            double phase = 0.0;
            advance(tc, phase, run_s);

            const Eigen::Vector2d before = pathPoint(tc.path(), phase);
            const double lap_before = tc.lapS();

            tc.setSizeMm(kMaxPathRadius * 1000.0);   // 20 -> 180 mm, floor moves
            const Eigen::Vector2d after = pathPoint(tc.path(), phase);

            // The floor really did move, or this test is not exercising the
            // indirect route it exists for.
            ASSERT_TRUE(tc.lapS() > lap_before);

            // Same ray from the centre, and the growth is the ratio asked for.
            ASSERT_TRUE(before.norm() > 1e-9 && after.norm() > 1e-9);
            ASSERT_NEAR(before.normalized().dot(after.normalized()), 1.0, 1e-12);
            ASSERT_NEAR(after.norm() / before.norm(), kMaxPathRadius / 0.02, 1e-9);
        }
    }
}

// And the lap the size slider leaves behind is the floor of the path it just
// made, not of the one it replaced.  **This is the stale-radius skew, asserted
// rather than argued**: the two lines inside `setSizeMm` are in one order and
// not the other, and reversing them fails here.
//
// It is measurable because the floor moves so far.  Every 20 mm shape sits at
// the two-second minimum; at 180 mm the floors are 4.52 s for the circle, 3.86
// for the square and 3.23 for the triangle.  Left at the small path's floor,
// the same shapes run their setpoint at 566, 482 and 403 mm/s against a 250
// mm/s cap — which is the bound that keeps the ball on the plate at all.
void test_growing_the_path_raises_the_lap_to_the_new_floor() {
    for (PathShape s : {PathShape::Circle, PathShape::Square, PathShape::Triangle}) {
        TrajectoryControls tc(shape(s, 0.02, 10.0, kSomeAccel));
        tc.setLapS(minPeriod(tc.path()));
        const double small_floor = tc.lapS();

        tc.setSizeMm(kMaxPathRadius * 1000.0);

        // The floor of the path it now HAS, to the last bit.
        ASSERT_NEAR(tc.lapS(), minPeriod(tc.path()), 1e-15);
        ASSERT_TRUE(tc.lapS() > small_floor + 1.0);

        // And what the other order would have left: the big path at the small
        // path's floor, which is 1.6 to 2.3 times over the cap the bound
        // exists to hold it under.
        ASSERT_TRUE(pathLength(tc.path()) / small_floor > 1.6 * kMaxSetpointSpeed);
        ASSERT_TRUE(pathLength(tc.path()) / tc.lapS() <= kMaxSetpointSpeed + 1e-12);
    }
}

// The same slider on a BLENDED path, which cannot be exactly radial and is not
// claimed to be: the fillet is set by speed rather than by size, so growing the
// path changes the fillet as well as the shape, and the setpoint slides very
// slightly around as well as outward.
//
// The promise the slider actually makes is that the target does not JUMP, and
// that survives: measured over both polygons at 200 phases, dragging 20 mm to
// 180 mm moves the setpoint off its own ray by at most 3.0 degrees and 0.8 mm.
// Against the 160 mm it travels outward, that is the target sliding to the new
// path rather than being thrown to a different part of it.
void test_the_size_slider_barely_rotates_a_filleted_setpoint() {
    for (PathShape s : {PathShape::Square, PathShape::Triangle}) {
        for (int i = 0; i < 200; ++i) {
            TrajectoryControls tc(shape(s, 0.02, 10.0, kSomeAccel));
            tc.setLapS(minPeriod(tc.path()));
            const double phase = i / 200.0;

            const Eigen::Vector2d before = pathPoint(tc.path(), phase);
            tc.setSizeMm(kMaxPathRadius * 1000.0);
            const Eigen::Vector2d after = pathPoint(tc.path(), phase);

            ASSERT_TRUE(before.norm() > 1e-9 && after.norm() > 1e-9);
            const double cos_between =
                before.normalized().dot(after.normalized());
            ASSERT_TRUE(cos_between > std::cos(3.5 * M_PI / 180.0));
            // And sideways, in millimetres, against a 20 mm path.
            const double sideways =
                (after.normalized() * before.norm() - before).norm();
            ASSERT_TRUE(sideways < 0.001);
        }
    }
}

// And the rate the phase advances at is the lap time, from the moment the
// slider moved — which is the whole of what the slider is entitled to do.
void test_a_lap_change_takes_effect_from_that_moment() {
    TrajectoryControls tc(shape(PathShape::Circle, 0.12, 10.0));
    double phase = 0.0;
    advance(tc, phase, 100.0);
    const double at_change = phase;

    tc.setLapS(5.0);
    advance(tc, phase, 1.0);
    // One second at a five-second lap is a fifth of a lap, not a tenth.
    double moved = phase - at_change;
    if (moved < 0.0) moved += 1.0;          // the lap wrapped in between
    ASSERT_NEAR(moved, 0.2, 1e-9);
}

// Changing shape must not move the setpoint either — the third control, and
// the one that was worst.  Phase is NOT comparable across shapes: the circle's
// phase zero is at +x while a polygon's first corner is at the top, so equal
// phase is a quarter of a lap apart.  Carrying it over threw the target 170 mm
// to the far side of a 120 mm path, and the loop hauled the ball across after
// it hard enough to drive the legs into the workspace clip.
void test_changing_the_shape_picks_up_where_the_setpoint_is() {
    const PathShape all[] = {PathShape::Circle, PathShape::Square,
                             PathShape::Triangle};
    double worst = 0.0;
    for (PathShape from : all) {
        for (PathShape to : all) {
            for (int i = 0; i < 64; ++i) {
                TrajectoryControls tc(shape(from, 0.12, 10.0));
                double phase = i / 64.0;
                frameLeftItAt(tc, phase);

                const Eigen::Vector2d before = tc.setpoint();
                tc.setShape(to, phase);
                const Eigen::Vector2d after = pathPoint(tc.path(), phase);
                worst = std::max(worst, (after - before).norm());

                // Switching to the SAME shape is exactly a no-op, at every
                // phase — the nearest point on a path to a point already on it
                // is that point.
                if (from == to) ASSERT_NEAR((after - before).norm(), 0.0, 1e-9);
            }
        }
    }
    // The shapes themselves differ, so a switch cannot be free — but it is
    // bounded by how far apart the two paths ARE, which is the least it could
    // possibly be.  Measured over every ordered pair at 64 phases, on a 120 mm
    // path: 35.15 mm circle <-> square, 42.43 mm square -> triangle, and 60.00
    // mm circle <-> triangle, which is exactly the r/2 a triangle's inradius
    // leaves.  Nothing like the 169.7 mm a carried phase gave at EVERY phase.
    ASSERT_TRUE(worst < 0.061);
}

// And the bound above is only meaningful next to what it replaced.  Carrying
// the phase across a shape change — the obvious simplification, and what the
// code did — is 170 mm on this path at EVERY phase in the lap.  Asserted so
// that reverting to it fails a test rather than merely feeling smoother.
void test_carrying_the_phase_across_a_shape_change_would_be_far_worse() {
    const SetpointPath circle = shape(PathShape::Circle, 0.12, 10.0);
    const SetpointPath square = shape(PathShape::Square, 0.12, 10.0);

    double worst_carried = 0.0, worst_reseeded = 0.0;
    for (int i = 0; i < 512; ++i) {
        const double u = i / 512.0;
        const Eigen::Vector2d on_circle = pathPoint(circle, u);
        // What carrying the phase does.
        worst_carried = std::max(worst_carried,
                                 (pathPoint(square, u) - on_circle).norm());
        // What re-seeding does.
        const Eigen::Vector2d seeded =
            pathPoint(square, phaseNearest(square, on_circle));
        worst_reseeded = std::max(worst_reseeded, (seeded - on_circle).norm());
    }
    ASSERT_TRUE(worst_carried > 0.16);      // measured 169.7 mm
    ASSERT_TRUE(worst_reseeded < 0.036);    // measured 35.1 mm
    ASSERT_TRUE(worst_reseeded * 4.0 < worst_carried);
}

// `phaseNearest` really is the nearest, not merely a good guess — it is solved
// in closed form, so it is worth checking against a brute-force sweep.
void test_the_nearest_phase_is_actually_the_nearest() {
    const Eigen::Vector2d probes[] = {
        Eigen::Vector2d(0.12, 0.0),   Eigen::Vector2d(-0.03, 0.09),
        Eigen::Vector2d(0.0, -0.2),   Eigen::Vector2d(0.001, 0.0),
        Eigen::Vector2d(0.4, 0.4),    Eigen::Vector2d(-0.07, -0.02),
    };
    // Blended as well as sharp: the filleted answer has to consider the arcs,
    // and the nearest point on an arc is very often in the middle of one — a
    // walk that only ever offered the straights would be wrong by up to the
    // fillet radius, which is 33 mm at the settings the sliders reach.
    for (double a : {0.0, kSomeAccel}) {
        for (PathShape s : {PathShape::Circle, PathShape::Square, PathShape::Triangle}) {
            const SetpointPath p = shape(s, 0.12, 4.0, a);
            for (const Eigen::Vector2d& t : probes) {
                const double got = (pathPoint(p, phaseNearest(p, t)) - t).norm();
                double brute = std::numeric_limits<double>::max();
                for (int i = 0; i < 20000; ++i)
                    brute = std::min(brute, (pathPoint(p, i / 20000.0) - t).norm());
                // The sweep can only ever be worse, up to its own resolution.
                ASSERT_TRUE(got <= brute + 1e-6);
            }
        }
    }
}

// The centre is equidistant from every point of a path, so it has no nearest
// one.  It answers phase zero rather than a NaN, which would poison the walk.
void test_the_centre_has_no_nearest_phase_and_says_so() {
    // Phase zero for every shape, and the same answer every time.  Left to the
    // polygon walk this came back 0.375 — the midpoint of whichever edge
    // rounding made shortest, which is a tie broken by noise.
    for (PathShape s : {PathShape::Circle, PathShape::Square, PathShape::Triangle}) {
        const double u = phaseNearest(shape(s, 0.12, 10.0), Eigen::Vector2d::Zero());
        ASSERT_TRUE(std::isfinite(u));
        ASSERT_NEAR(u, 0.0, 1e-15);
    }
    // And a `Fixed` path is a point, not a lap: every phase is that point.
    ASSERT_NEAR(phaseNearest(shape(PathShape::Fixed), Eigen::Vector2d(0.1, 0.1)),
                0.0, 1e-15);
}

// Drawing: a sharp polygon is drawn by its corners, a circle by samples.  A
// square drawn as a chord of samples is a square drawn wrong.
void test_the_outline_is_drawable() {
    Eigen::Matrix<double, 2, Eigen::Dynamic> pts;
    pathOutline(shape(PathShape::Square), 64, &pts);
    ASSERT_EQ((int)pts.cols(), 5);                     // four corners, closed
    ASSERT_NEAR((pts.col(0) - pts.col(4)).norm(), 0.0, 1e-12);

    pathOutline(shape(PathShape::Triangle), 64, &pts);
    ASSERT_EQ((int)pts.cols(), 4);

    pathOutline(shape(PathShape::Circle), 64, &pts);
    ASSERT_EQ((int)pts.cols(), 65);
    ASSERT_NEAR(pts.col(0).norm(), 0.12, 1e-12);

    pathOutline(shape(PathShape::Fixed), 64, &pts);
    ASSERT_EQ((int)pts.cols(), 0);
}

// And a filleted one is drawn filleted.  A square outlined with sharp corners
// over a setpoint that rounds them is a picture of a path the ball is not
// being sent round — the visitor would read the gap at the corner as the
// controller failing rather than as the corner not being there.
//
// The claim is not "there are more points"; it is that every point drawn is ON
// the path.  Checked against `phaseNearest`, which is the path's own answer to
// "where is the closest bit of me".
void test_the_outline_draws_the_path_the_setpoint_runs() {
    Eigen::Matrix<double, 2, Eigen::Dynamic> pts;
    for (PathShape s : {PathShape::Square, PathShape::Triangle}) {
        const SetpointPath p = shape(s, kMaxPathRadius, 4.0, kSomeAccel);
        pathOutline(p, 64, &pts);
        ASSERT_TRUE(pts.cols() > pathCorners(s) + 1);        // fillets sampled
        ASSERT_NEAR((pts.col(0) - pts.col(pts.cols() - 1)).norm(), 0.0, 1e-12);
        for (int i = 0; i < pts.cols(); ++i) {
            const Eigen::Vector2d at = pts.col(i);
            const Eigen::Vector2d on = pathPoint(p, phaseNearest(p, at));
            ASSERT_NEAR((at - on).norm(), 0.0, 1e-9);
        }
        // Drawn with a sharp outline instead, the corners would stand up to
        // 10 mm proud of the path the setpoint runs.  That is the picture this
        // test exists to forbid.
        SetpointPath sharp = p;
        sharp.accel_max = 0.0;
        Eigen::Matrix<double, 2, Eigen::Dynamic> sharp_pts;
        pathOutline(sharp, 64, &sharp_pts);
        double worst = 0.0;
        for (int i = 0; i < sharp_pts.cols(); ++i) {
            const Eigen::Vector2d at = sharp_pts.col(i);
            worst = std::max(worst, (at - pathPoint(p, phaseNearest(p, at))).norm());
        }
        ASSERT_TRUE(worst > 0.005);
    }
}

// ---------------------------------------------------------------------------
// The controls, through the application's own frame (#28)
// ---------------------------------------------------------------------------
//
// Everything above drives `TrajectoryControls` and reads the path back.  That
// is the panel's object, but it is not yet the panel's FRAME: the phase lives
// in `SimState` and `stepSim` advances it, so a test that never runs a step
// cannot see the two wired to each other.
//
// **That gap is what this ticket is about.**  Both of the last two defects were
// a well-tested rule called from a place nothing tested — the review of #24
// round two said in as many words that `plate_view.cpp` could go back to
// deriving the phase from the clock and every test would still pass, and round
// three then shipped a shape change that threw the target 169.7 mm because the
// combo handler was view code.  So these run frames.
//
// The plate is the shipped one, which means `accel_max` is the real 1.888
// m/s^2 rather than a number this file chose.  The setpoint tests above cannot
// use it — a filleted path changes shape when the lap changes, so "the setpoint
// does not move" stops being exact — and the circle is the way through: it has
// no corners, so it is the one shape whose geometry the fillet does not touch.

/// The shipped cascade's parameters — the preset the application opens on.
const std::vector<PhysicalParam>& cascadeParams() {
    static const std::vector<ModelEntry> models = getBuiltinModels();
    static const std::vector<PhysicalParam> params = cascadeModel(models).params;
    return params;
}

/// The three objects `PlateView` holds, and the two calls it makes either side
/// of `stepSim`.
///
/// It carries no rule of its own — not the read-then-advance order, not the lap
/// floor, not the re-seed.  Those are `stepSim`'s and `TrajectoryControls`'.
/// What it is for is that they are wired to one another, which is the one thing
/// neither a test of the path arithmetic nor a test of the controls can see.
///
/// No gain and no ball.  What is being measured is where the setpoint went, and
/// nothing in the setpoint stage depends on who is driving the legs; a ball
/// nobody simulates cannot roll off partway through and leave the rest of a
/// hundred-second run measuring an empty plate.  `test_trajectory` is where the
/// ball's side of this is pinned.
struct Panel {
    SimPlate plate;
    TrajectoryControls traj;
    SimState sim;
    SimInput in;

    explicit Panel(const SetpointPath& opening)
        : plate(cascadePlate(cascadeParams())),
          traj(plate.feasible(opening)),
          sim(simStart(plate, cascadeDesign(cascadeParams()).home_leg_rad)) {
        in.design = cascadeDesign(cascadeParams());
        in.closed_loop = false;
        in.ball_enabled = false;
    }

    /// One frame, exactly as `PlateView::step` runs it.
    Eigen::Vector2d step() {
        traj.toInput(in);
        const SimReport r = stepSim(plate, in, sim);
        traj.fromReport(r);
        return r.setpoint;
    }

    void run(double seconds) {
        for (int i = 0; i < static_cast<int>(seconds / in.dt); ++i) step();
    }
};

// The frame drives the path the panel is holding, and hands the answer back to
// it.  Both halves, because each can be broken on its own: a panel whose path
// is not the one being run draws an outline the ball is not chasing, and a
// panel that does not read the setpoint back shows the visitor a stale target
// and re-seeds the next shape change from it.
void test_a_frame_drives_the_panel_s_own_path_and_hands_it_back() {
    Panel p(shape(PathShape::Circle, 0.12, 10.0));

    // The fillet the step will stamp is the fillet the panel already holds, so
    // the outline it draws and the lap floor it offers are about the path that
    // is actually being run.
    ASSERT_NEAR(p.traj.path().accel_max, p.plate.maxBallAccel(), 1e-15);

    const double phase_before = p.sim.path_phase;
    const Eigen::Vector2d was = p.traj.setpoint();
    const Eigen::Vector2d drove = p.step();

    // Where the step drove the setpoint is where the panel's own path says it
    // should be, at the phase the frame opened on — read then advance.
    ASSERT_NEAR((drove - pathPoint(p.traj.path(), phase_before)).norm(), 0.0, 1e-15);
    ASSERT_TRUE(p.sim.path_phase > phase_before);
    // And the panel is holding it.  It started at the centre, so this also says
    // the frame moved it at all.
    ASSERT_NEAR((p.traj.setpoint() - drove).norm(), 0.0, 1e-15);
    ASSERT_TRUE((p.traj.setpoint() - was).norm() > 0.1);
}

// A held setpoint is the visitor's, and a frame does not touch it — the other
// half of the same wiring.  Drag it, run a minute, find it where it was left.
void test_a_held_setpoint_is_the_visitor_s_and_a_frame_leaves_it_alone() {
    Panel p{SetpointPath{}};                 // PathShape::Fixed
    p.traj.setHeldSetpointMm(40.0, -25.0);

    const Eigen::Vector2d held = p.traj.setpoint();
    for (int i = 0; i < 3600; ++i)
        ASSERT_NEAR((p.step() - held).norm(), 0.0, 1e-15);
    ASSERT_NEAR((p.traj.setpoint() - held).norm(), 0.0, 1e-15);
    // And the lap it is not running has not moved under it, so choosing a shape
    // resumes where the phase was left rather than wherever the clock got to.
    ASSERT_NEAR(p.sim.path_phase, 0.0, 1e-15);
}

// **Criterion: deriving the phase from absolute time again fails a test.**
//
// The phase used to be `t / period_s`.  Changing the lap moves that quantity by
// `t dT / T^2`, and `t` is the whole time the demo has been running — measured
// at 100 s, a nudge from 10.0 to 9.5 s threw the setpoint 170 degrees round the
// path.  `test_changing_the_lap_leaves_the_setpoint_where_it_is` says the rule;
// this says the application obeys it, which is the part that was untested.
//
// **Where the revert would now have to be made, stated so the attribution is
// not overclaimed.**  Since #30 the phase is `stepSim`'s, so what goes red here
// is the step — not this ticket's seam.  The panel-side half of the original
// defect is gone by construction rather than by test: `PlateView` holds no
// phase and reads no clock for one, so there is nothing left in it to derive a
// phase from.  The criterion is kept because the property is the one the review
// of #24 round two said was reachable from `plate_view.cpp` with every test
// still passing; what has changed is which file it is reachable from.
void test_the_lap_slider_does_not_move_the_setpoint_through_the_frame() {
    Panel p(shape(PathShape::Circle, 0.12, 10.0));
    p.run(100.0);

    const Eigen::Vector2d before = p.traj.setpoint();
    p.traj.setLapS(9.5);
    const Eigen::Vector2d after = p.step();

    // One frame's worth of travel and no more.  At 75 mm/s that is 1.3 mm; the
    // defect was 238 mm on this path, so nothing here turns on the bound.
    ASSERT_TRUE((after - before).norm() < 0.002);

    // Said the other way, so that the margin is not the whole of the claim:
    // the setpoint is nowhere near where the clock would have put it.
    const Eigen::Vector2d by_the_clock = pathPoint(p.traj.path(), 100.0 / 9.5);
    ASSERT_TRUE((after - by_the_clock).norm() > 0.1);
}

// **Criterion: carrying the phase across a shape change again fails a test.**
//
// The instance that reached a user.  The combo is view code, so nothing could
// see it; `setShape` is not, and this drives it through the same frame the
// panel does.
void test_a_shape_change_through_the_frame_does_not_throw_the_setpoint() {
    for (PathShape to : {PathShape::Square, PathShape::Triangle, PathShape::Circle}) {
        Panel p(shape(PathShape::Circle, 0.12, 10.0));
        p.run(37.0);                          // no particular phase

        const Eigen::Vector2d before = p.traj.setpoint();
        const double carried = p.sim.path_phase;
        p.traj.setShape(to, p.sim.path_phase);
        const Eigen::Vector2d after = p.step();

        // At most the two shapes are apart, plus the frame's own travel: 35 mm
        // to the square, 60 mm to the triangle, and nothing at all back to the
        // circle it is already on.
        ASSERT_TRUE((after - before).norm() < 0.062);
        if (to == PathShape::Circle)
            ASSERT_TRUE((after - before).norm() < 0.002);

        // And what carrying the phase would have done instead, on the path the
        // step is now running: the far side of it.  At the phase 37 seconds
        // leaves, re-seeding moves the target 23 mm to the square and 54 mm to
        // the triangle, against 161 mm and 150 mm carried.
        if (to != PathShape::Circle) {
            const Eigen::Vector2d kept = pathPoint(p.traj.path(), carried);
            ASSERT_TRUE((kept - before).norm() > 0.14);
            ASSERT_TRUE((kept - before).norm() > 2.5 * (after - before).norm());
        }
    }
}

// **Criterion: applying the lap floor against a stale radius again fails a
// test** — through the frame, where the consequence is a speed rather than a
// number.  The size slider is dragged to the maximum on a path already at its
// floor, and what comes out has to be a setpoint the cap permits, measured by
// how far it actually travels in a frame rather than by asking the path.
void test_growing_the_path_mid_run_keeps_the_setpoint_under_the_cap() {
    Panel p(shape(PathShape::Circle, 0.02, 10.0));
    p.traj.setLapS(minPeriod(p.traj.path()));
    p.run(5.0);

    p.traj.setSizeMm(kMaxPathRadius * 1000.0);

    // A lap of the new path, checked frame by frame.  The first frame after the
    // drag is skipped: the setpoint slides outward to the bigger path on it,
    // which is the size slider doing its job and not the speed of the lap.
    Eigen::Vector2d last = p.step();
    double fastest = 0.0;
    for (int i = 0; i < static_cast<int>(p.traj.lapS() / p.in.dt); ++i) {
        const Eigen::Vector2d now = p.step();
        fastest = std::max(fastest, (now - last).norm() / p.in.dt);
        last = now;
    }
    ASSERT_TRUE(fastest <= kMaxSetpointSpeed + 1e-6);
    // And it is running near the cap rather than crawling, or the bound would
    // not be the thing this test is about.
    ASSERT_TRUE(fastest > 0.9 * kMaxSetpointSpeed);
}

// ---------------------------------------------------------------------------
// Where `accel_max` comes from (#31)
// ---------------------------------------------------------------------------
//
// Everything above makes up its own number for what the plate can give the
// ball.  These do not: the fillet is only feasible if the number is the real
// one, and "derived from the plant rather than chosen" is the claim that makes
// it a bound instead of a taste.

// `a_max = (5/7) g sin(theta_max)`, and `theta_max` is swept rather than
// written down.  The three parts are checked separately because each is a
// different way of being wrong.
void test_a_max_is_the_ball_s_own_acceleration_at_the_swept_tilt() {
    const auto models = getBuiltinModels();
    const auto& e = cascadeModel(models);
    const TableParams tp = cascadeMechanism(e.params);
    const double g = cascadeGravity(e.params);
    const double home = cascadeHomeLegAngle(e.params);
    const TableKinematics tk(tp);

    const double theta_max =
        tk.max_conditioned_tilt(tk.home_pose(home).z_c, kRatesUntrustworthyAbove);
    // Measured on the shipped geometry: 15.63 degrees, 1.888 m/s^2.
    ASSERT_NEAR(theta_max * 180.0 / M_PI, 15.63, 0.05);
    ASSERT_NEAR(maxBallAccel(tk, g, home),
                RollingBallDynamics::rolling_factor() * g * std::sin(theta_max),
                1e-15);
    ASSERT_NEAR(maxBallAccel(tk, g, home), 1.888, 0.01);

    // And it is the plate's, not a constant: `cascadePlate` measures the same
    // thing the same way, which is what lets `stepSim` stamp it.
    ASSERT_NEAR(cascadePlate(e.params).maxBallAccel(), maxBallAccel(tk, g, home),
                1e-15);
}

// The tilt really is the largest one the plate can hold and be believed at:
// the predicate holds there and fails just past it.  Both halves, because a
// sweep that returned zero would pass either one alone.
void test_the_swept_tilt_is_the_last_one_that_holds() {
    const auto models = getBuiltinModels();
    const auto& e = cascadeModel(models);
    const TableKinematics tk(cascadeMechanism(e.params));
    const double z_c = tk.home_pose(cascadeHomeLegAngle(e.params)).z_c;
    const double theta_max =
        tk.max_conditioned_tilt(z_c, kRatesUntrustworthyAbove);
    ASSERT_TRUE(theta_max > 0.0);

    // The two halves of the predicate, asked separately, so that which of them
    // binds is something this file can state rather than assume.
    const auto reachable = [&](double tilt) {
        for (int i = 0; i < 180; ++i) {
            const TablePose pose =
                TableKinematics::tilted_pose(tilt, 2.0 * M_PI * i / 180.0, z_c);
            for (int leg = 0; leg < 3; ++leg) {
                const auto a = tk.inverse_kinematics_leg(leg, pose);
                if (!a) return false;
                if (*a < tk.params().alpha_min || *a > tk.params().alpha_max)
                    return false;
            }
        }
        return true;
    };
    const auto trusted = [&](double tilt) {
        for (int i = 0; i < 180; ++i) {
            const TablePose pose =
                TableKinematics::tilted_pose(tilt, 2.0 * M_PI * i / 180.0, z_c);
            const IKResult ik = tk.inverse_kinematics(pose);
            if (!(tk.condition_number(ik.alpha, pose) < kRatesUntrustworthyAbove))
                return false;
        }
        return true;
    };

    ASSERT_TRUE(reachable(theta_max) && trusted(theta_max));
    // A tenth of a degree past it, one of the two has given way.
    const double past = theta_max + 0.1 * M_PI / 180.0;
    ASSERT_TRUE(!reachable(past) || !trusted(past));

    // And the threshold is the one this repository already argued for, not a
    // second opinion about when the mechanism is in trouble.
    ASSERT_NEAR(kRatesUntrustworthyAbove, 20.0, 1e-15);
}

// **Which of the two binds is a finding, and #31 expected the other one.**
//
// The ticket says "the workspace maximum is explicitly not the answer: the
// workspace edge is precisely where the condition number blows up".  Measured,
// the blow-up is real and it is 1.5 degrees OUTSIDE the reachable set, so on
// the shipped plate the travel binds first and `a_max` IS the workspace
// maximum.  Asserted rather than written in a comment, because it is the sort
// of claim that quietly stops being true.
void test_on_the_shipped_plate_it_is_the_travel_that_binds() {
    const auto models = getBuiltinModels();
    const auto& e = cascadeModel(models);
    const double home = cascadeHomeLegAngle(e.params);
    const TableKinematics tk(cascadeMechanism(e.params));
    const double z_c = tk.home_pose(home).z_c;

    // Where each criterion would stop, on its own.
    const double reach = tk.max_conditioned_tilt(z_c, 1e30);
    const double both = tk.max_conditioned_tilt(z_c, kRatesUntrustworthyAbove);
    // Measured: 15.63 degrees either way — the condition line is never
    // reached, so raising the limit to absurdity changes nothing.
    ASSERT_NEAR(both, reach, 1e-12);
    ASSERT_NEAR(both * 180.0 / M_PI, 15.63, 0.05);
}

// Which is exactly why the gate has to be tested somewhere it DOES bind, or
// this suite would go on passing with the condition test deleted.
//
// Two places it binds, both measured.  A lower limit bites at once on the
// shipped plate — the criterion is live, it is simply not the tighter of the
// two at 20.  And at 210 mm legs it binds outright at the shipped threshold:
// the plate reaches 21.80 degrees and is only to be believed to 20.63.  That
// is the property the gate exists for — a longer leg must not buy acceleration
// by reaching into rates nobody should trust.
void test_the_condition_gate_binds_where_the_mechanism_is_weaker() {
    const auto models = getBuiltinModels();
    const auto& e = cascadeModel(models);
    const double home = cascadeHomeLegAngle(e.params);

    const TableKinematics shipped(cascadeMechanism(e.params));
    const double z_c = shipped.home_pose(home).z_c;
    const double at_20 = shipped.max_conditioned_tilt(z_c, kRatesUntrustworthyAbove);
    // Measured: 15.51 degrees at a limit of 15, 13.66 at 10.
    ASSERT_TRUE(shipped.max_conditioned_tilt(z_c, 15.0) < at_20);
    ASSERT_TRUE(shipped.max_conditioned_tilt(z_c, 10.0) <
                shipped.max_conditioned_tilt(z_c, 15.0));

    TableParams longer = cascadeMechanism(e.params);
    longer.L1 = 0.21;
    longer.L2 = 0.21;
    const TableKinematics big(longer);
    const double big_z = big.home_pose(home).z_c;
    const double trusted = big.max_conditioned_tilt(big_z, kRatesUntrustworthyAbove);
    const double reach = big.max_conditioned_tilt(big_z, 1e30);
    // The gate really is the tighter of the two here, by over a degree.
    ASSERT_TRUE(trusted < reach - 0.5 * M_PI / 180.0);
    ASSERT_NEAR(trusted * 180.0 / M_PI, 20.68, 0.1);
    ASSERT_NEAR(reach * 180.0 / M_PI, 21.80, 0.1);
}

// A derived property of the plant moves when the plant does — which is the
// behaviour every other bound in `setpoint_path.h` already has, and the reason
// `theta_max` is swept rather than written as a literal.  A longer leg reaches
// further over, so the ball can be accelerated harder, so the corners need
// less rounding.
void test_a_max_moves_with_the_legs_and_the_fillet_moves_with_it() {
    const auto models = getBuiltinModels();
    const auto& e = cascadeModel(models);
    const double g = cascadeGravity(e.params);
    const double home = cascadeHomeLegAngle(e.params);

    double previous_accel = 0.0, previous_fillet = 1.0;
    // Measured: 1.52, 1.89, 2.25, 2.47 m/s^2.
    for (double leg : {0.12, 0.15, 0.18, 0.21}) {
        TableParams tp = cascadeMechanism(e.params);
        tp.L1 = leg;
        tp.L2 = leg;
        const double a = maxBallAccel(TableKinematics(tp), g, home);
        ASSERT_TRUE(a > previous_accel);
        previous_accel = a;

        const SetpointPath p = shape(PathShape::Square, kMaxPathRadius, 4.0, a);
        ASSERT_TRUE(filletRadius(p) < previous_fillet);
        previous_fillet = filletRadius(p);
    }
    // The whole range is a real spread rather than a rounding difference.
    ASSERT_TRUE(previous_accel > 1.5 * 1.52);
}

}  // namespace

int main() {
    test_every_path_closes();
    test_no_path_leaves_its_radius();
    test_every_shape_reaches_its_radius();
    test_a_filleted_polygon_falls_short_of_its_radius_by_the_fillet();
    test_the_polygons_are_walked_at_constant_speed();
    test_velocity_is_the_derivative_of_position();
    test_the_fillet_keeps_the_velocity_the_position_s_derivative_at_the_corners();
    test_a_sharp_corner_is_a_step_in_the_reference_velocity();
    test_a_filleted_corner_turns_the_reference_velocity_continuously();
    test_the_fillet_radius_is_v_squared_over_a_max();
    test_the_fillet_grows_as_the_lap_tightens();
    test_the_filleted_reference_never_exceeds_a_max();
    test_the_fillet_is_capped_at_the_incircle();
    test_a_filleted_polygon_is_shorter_and_so_is_its_fastest_lap();
    test_the_degenerate_cases_behave();
    test_changing_the_lap_leaves_the_setpoint_where_it_is();
    test_the_reference_velocity_steps_on_a_lap_change();
    test_changing_the_size_moves_the_setpoint_radially_only();
    test_growing_the_path_raises_the_lap_to_the_new_floor();
    test_the_size_slider_barely_rotates_a_filleted_setpoint();
    test_a_lap_change_takes_effect_from_that_moment();
    test_changing_the_shape_picks_up_where_the_setpoint_is();
    test_carrying_the_phase_across_a_shape_change_would_be_far_worse();
    test_the_nearest_phase_is_actually_the_nearest();
    test_the_centre_has_no_nearest_phase_and_says_so();
    test_the_outline_is_drawable();
    test_the_outline_draws_the_path_the_setpoint_runs();
    test_a_frame_drives_the_panel_s_own_path_and_hands_it_back();
    test_a_held_setpoint_is_the_visitor_s_and_a_frame_leaves_it_alone();
    test_the_lap_slider_does_not_move_the_setpoint_through_the_frame();
    test_a_shape_change_through_the_frame_does_not_throw_the_setpoint();
    test_growing_the_path_mid_run_keeps_the_setpoint_under_the_cap();
    test_a_max_is_the_ball_s_own_acceleration_at_the_swept_tilt();
    test_the_swept_tilt_is_the_last_one_that_holds();
    test_on_the_shipped_plate_it_is_the_travel_that_binds();
    test_the_condition_gate_binds_where_the_mechanism_is_weaker();
    test_a_max_moves_with_the_legs_and_the_fillet_moves_with_it();
    std::printf("test_setpoint_path: all passed\n");
    return 0;
}
