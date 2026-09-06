// tests/test_setpoint_path.cpp
//
// The paths the ball is asked to trace.  Nothing here touches the controller —
// these are claims about the shapes themselves, and they are worth making
// separately because a tracking failure is much easier to read when the thing
// being tracked is known to be right.
#include "setpoint_path.h"

#include "test_helpers.h"

#include <algorithm>
#include <cmath>
#include <limits>

using namespace caliburn;

namespace {

SetpointPath shape(PathShape s, double r = 0.12, double period = 10.0) {
    SetpointPath p;
    p.shape = s;
    p.radius_m = r;
    p.period_s = period;
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

// Constant speed, not constant angle.  A polygon traversed by sweeping an
// angle would crawl through its corners, which is precisely where the
// interesting behaviour is.
void test_the_polygons_are_walked_at_constant_speed() {
    for (PathShape s : {PathShape::Square, PathShape::Triangle}) {
        const SetpointPath p = shape(s, 0.12, 6.0);
        const double expect = pathLength(p) / p.period_s;
        for (int i = 0; i < 600; ++i) {
            const double u = (i + 0.5) / 600.0;         // off the corners
            ASSERT_NEAR(pathVelocity(p, u).norm(), expect, 1e-9);
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

// A corner IS a step in the reference velocity — that is what a corner means,
// and it is why the ball rounds one.  Worth asserting rather than assuming:
// if the path ever stopped having corners the demo would lose its point.
void test_a_corner_is_a_step_in_the_reference_velocity() {
    const SetpointPath p = shape(PathShape::Square, 0.12, 8.0);
    const double corner_u = 0.25;                  // end of the first edge
    const Eigen::Vector2d before = pathVelocity(p, corner_u - 1e-5);
    const Eigen::Vector2d after = pathVelocity(p, corner_u + 1e-5);
    // A square's corner turns the velocity through a right angle.
    ASSERT_NEAR(before.normalized().dot(after.normalized()), 0.0, 1e-6);
    ASSERT_NEAR(before.norm(), after.norm(), 1e-9);
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

/// `plate_view`'s own arrangement: a path, an accumulated phase, and the two
/// sliders applied the way the panel applies them.
struct Run {
    SetpointPath path;
    double phase = 0.0;

    void advance(double seconds) {
        const double dt = 1.0 / 60.0;
        for (int i = 0; i < static_cast<int>(seconds / dt); ++i)
            phase = stepPath(path, phase, dt).next_phase;
    }

    Eigen::Vector2d point() const { return pathPoint(path, phase); }
    Eigen::Vector2d velocity() const { return pathVelocity(path, phase); }

    /// The lap slider.  Held off the floor by `clampPeriod`, which is the
    /// panel's own rule rather than a copy of it — see `setpoint_path.h`.
    void setLap(double period_s) { path.period_s = clampPeriod(path, period_s); }

    /// The size slider — which moves the lap floor with it, and so was the
    /// second way into the same bug: raising the radius raises the floor,
    /// which pushes the lap up, which was a change of `period_s`.
    void setSize(double radius_m) {
        path.radius_m = radius_m;
        path.period_s = clampPeriod(path, path.period_s);
    }

    /// The shape combo.  The new shape picks up nearest to where the setpoint
    /// already is — `phaseNearest`, the panel's own rule — because phase means
    /// a different place on each shape.
    void setShape(PathShape s) {
        const Eigen::Vector2d was = point();
        path.shape = s;
        path.period_s = clampPeriod(path, path.period_s);
        phase = phaseNearest(path, was);
    }
};

// The lap slider changes the RATE and nothing else.  Not approximately —
// exactly: the phase is a number the slider does not touch.
void test_changing_the_lap_leaves_the_setpoint_where_it_is() {
    for (PathShape s : {PathShape::Circle, PathShape::Square, PathShape::Triangle}) {
        // Deliberately at four run lengths spanning four orders of magnitude.
        // The old failure was proportional to `t`, so a test that only ever
        // ran for a second would have passed against the bug.
        for (double run_s : {0.5, 10.0, 100.0, 1000.0}) {
            Run r;
            r.path = shape(s, 0.12, 10.0);
            r.advance(run_s);

            const Eigen::Vector2d before = r.point();
            r.setLap(9.5);
            const Eigen::Vector2d after = r.point();
            ASSERT_NEAR((after - before).norm(), 0.0, 1e-15);

            // And again the other way, and by a lot rather than a nudge.
            const Eigen::Vector2d before_2 = r.point();
            r.setLap(30.0);
            ASSERT_NEAR((r.point() - before_2).norm(), 0.0, 1e-15);
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
        Run r;
        r.path = shape(s, 0.12, 10.0);
        r.advance(100.0);

        const Eigen::Vector2d before = r.velocity();
        r.setLap(5.0);                      // above every shape's floor at 120 mm
        const Eigen::Vector2d after = r.velocity();
        ASSERT_NEAR(r.path.period_s, 5.0, 1e-12);

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
            Run r;
            r.path = shape(s, 0.02, 10.0);
            r.setLap(minPeriod(r.path));          // sitting on the floor
            r.advance(run_s);

            const Eigen::Vector2d before = r.point();
            const double lap_before = r.path.period_s;

            r.setSize(kMaxPathRadius);            // 20 mm -> 180 mm, floor moves
            const Eigen::Vector2d after = r.point();

            // The floor really did move, or this test is not exercising the
            // indirect route it exists for.
            ASSERT_TRUE(r.path.period_s > lap_before);

            // Same ray from the centre, and the growth is the ratio asked for.
            ASSERT_TRUE(before.norm() > 1e-9 && after.norm() > 1e-9);
            ASSERT_NEAR(before.normalized().dot(after.normalized()), 1.0, 1e-12);
            ASSERT_NEAR(after.norm() / before.norm(), kMaxPathRadius / 0.02, 1e-9);
        }
    }
}

// And the rate the phase advances at is the lap time, from the moment the
// slider moved — which is the whole of what the slider is entitled to do.
void test_a_lap_change_takes_effect_from_that_moment() {
    Run r;
    r.path = shape(PathShape::Circle, 0.12, 10.0);
    r.advance(100.0);
    const double at_change = r.phase;

    r.setLap(5.0);
    r.advance(1.0);
    // One second at a five-second lap is a fifth of a lap, not a tenth.
    double moved = r.phase - at_change;
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
                Run r;
                r.path = shape(from, 0.12, 10.0);
                r.phase = i / 64.0;

                const Eigen::Vector2d before = r.point();
                r.setShape(to);
                const Eigen::Vector2d after = r.point();
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
    for (PathShape s : {PathShape::Circle, PathShape::Square, PathShape::Triangle}) {
        const SetpointPath p = shape(s, 0.12, 10.0);
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

// Drawing: a polygon is drawn by its corners, a circle by samples.  A square
// drawn as a chord of samples is a square drawn wrong.
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

}  // namespace

int main() {
    test_every_path_closes();
    test_no_path_leaves_its_radius();
    test_every_shape_reaches_its_radius();
    test_the_polygons_are_walked_at_constant_speed();
    test_velocity_is_the_derivative_of_position();
    test_a_corner_is_a_step_in_the_reference_velocity();
    test_the_degenerate_cases_behave();
    test_changing_the_lap_leaves_the_setpoint_where_it_is();
    test_the_reference_velocity_steps_on_a_lap_change();
    test_changing_the_size_moves_the_setpoint_radially_only();
    test_a_lap_change_takes_effect_from_that_moment();
    test_changing_the_shape_picks_up_where_the_setpoint_is();
    test_carrying_the_phase_across_a_shape_change_would_be_far_worse();
    test_the_nearest_phase_is_actually_the_nearest();
    test_the_centre_has_no_nearest_phase_and_says_so();
    test_the_outline_is_drawable();
    std::printf("test_setpoint_path: all passed\n");
    return 0;
}
