// tests/test_ball_contact.cpp
//
// Whether a ball stays on the plate is a question about the surface's
// ACCELERATION, and the rolling model cannot ask it: `RollingBallDynamics`
// carries [x, y, vx, vy] in the plate frame, has no normal force, and glues
// the ball down forever.  The plate heaves hard — `z_c = 0.30 sin(alpha)` for
// this geometry, so a 20 degree leg swing moves the table 74 mm vertically and
// the loop does that in about a tenth of a second — so the glue is doing real
// work.
//
// These tests build the normal force up from cases where the answer is known
// by hand, then point it at the shipped tuning.  See issue #23.
#include "ball_contact.h"

#include "cascade_fixture.h"
#include "test_helpers.h"

#include <cmath>

using namespace caliburn;

namespace {

constexpr double kDeg = M_PI / 180.0;
constexpr double kG = 9.81;
constexpr double kR = kFixtureBallRadius;

TableParams plate() {
    return cascadeMechanism(cascadeModel(getBuiltinModels()).params);
}

/// A plate that is not moving at all, at a given tilt.
PlateMotion still(const TableKinematics& tk, double phi, double theta, double z_c) {
    PlateMotion m;
    m.R = tk.table_rotation(phi, theta);
    m.c = Eigen::Vector3d(0, 0, z_c);
    return m;
}

// A ball sitting on a level, motionless plate is held by exactly its weight.
// The whole calculation reduces to this, and if it does not, nothing further
// down is worth reading.
void test_a_still_level_plate_holds_the_ball_at_one_g() {
    const TableKinematics tk(plate());
    const PlateMotion m = still(tk, 0.0, 0.0, 0.2121);
    const double N = normalAccel(m, Eigen::Vector3d(0.0, 0.0, kR),
                                 Eigen::Vector2d::Zero(), kG);
    ASSERT_NEAR(N, kG, 1e-12);
}

// Tilt it and the normal takes only gravity's share along the normal, which is
// the cosine.  The rest is what rolls the ball, and belongs to the roller.
void test_a_tilted_still_plate_holds_the_cosine_share() {
    const TableKinematics tk(plate());
    for (double deg : {5.0, 15.0, 30.0}) {
        const PlateMotion m = still(tk, deg * kDeg, 0.0, 0.2121);
        const double N = normalAccel(m, Eigen::Vector3d(0.0, 0.0, kR),
                                     Eigen::Vector2d::Zero(), kG);
        ASSERT_NEAR(N, kG * std::cos(deg * kDeg), 1e-12);
    }
}

// The case the whole ticket is about: drop the plate away from under the ball.
// At exactly one g of downward heave the ball is weightless — the plate is
// falling with it and touching it with nothing.  Below that, contact is over.
void test_heave_cancels_the_normal_force_at_one_g() {
    const TableKinematics tk(plate());

    auto heaving = [&](double zc_ddot) {
        PlateMotion m = still(tk, 0, 0, 0.2121);
        m.c_ddot = Eigen::Vector3d(0, 0, zc_ddot);
        return normalAccel(m, Eigen::Vector3d(0, 0, kR),
                           Eigen::Vector2d::Zero(), kG);
    };

    ASSERT_NEAR(heaving(0.0), kG, 1e-12);
    ASSERT_NEAR(heaving(-kG), 0.0, 1e-9);        // weightless, exactly
    ASSERT_TRUE(heaving(-1.5 * kG) < 0.0);       // separated
    ASSERT_NEAR(heaving(kG), 2.0 * kG, 1e-9);    // and pushed twice as hard
}

// Steady rotation is centripetal only, and about a horizontal axis through the
// table centre that pulls straight down on the ball's own radius — the same
// for a ball anywhere on the plate, because a sideways offset adds nothing to
// the vertical component.  Worth pinning precisely: it is the term most easily
// confused with the next one.
void test_steady_rotation_pulls_on_the_balls_own_radius() {
    const TableKinematics tk(plate());

    PlateMotion m = still(tk, 0, 0, 0.2121);
    m.omega = Eigen::Vector3d(2.0, 0.0, 0.0);   // steady: omega_dot stays zero

    const double centred = normalAccel(m, Eigen::Vector3d(0, 0, kR),
                                       Eigen::Vector2d::Zero(), kG);
    const double off = normalAccel(m, Eigen::Vector3d(0, 0.10, kR),
                                   Eigen::Vector2d::Zero(), kG);
    ASSERT_NEAR(centred, kG - 2.0 * 2.0 * kR, 1e-9);
    ASSERT_NEAR(off, centred, 1e-12);
}

// Angular ACCELERATION is the one that cares where the ball is: the plate
// tipping about its centre lifts one side and drops the other, so two balls on
// opposite sides are held by different forces at the same instant.  This is the
// term that takes the ball off the plate when the loop slams the legs over.
void test_angular_acceleration_lifts_one_side_and_drops_the_other() {
    const TableKinematics tk(plate());
    const double alpha_dd = 40.0;   // [rad/s^2] about x

    PlateMotion m = still(tk, 0, 0, 0.2121);
    m.omega_dot = Eigen::Vector3d(alpha_dd, 0.0, 0.0);   // omega itself is zero

    const double centred = normalAccel(m, Eigen::Vector3d(0, 0, kR),
                                       Eigen::Vector2d::Zero(), kG);
    const double plus_y = normalAccel(m, Eigen::Vector3d(0, 0.10, kR),
                                      Eigen::Vector2d::Zero(), kG);
    const double minus_y = normalAccel(m, Eigen::Vector3d(0, -0.10, kR),
                                       Eigen::Vector2d::Zero(), kG);

    // omega_dot x r has vertical component alpha_dd * y, so the two sides
    // straddle the centred case by that much, symmetrically.
    ASSERT_NEAR(plus_y - centred, alpha_dd * 0.10, 1e-6);
    ASSERT_NEAR(centred - minus_y, alpha_dd * 0.10, 1e-6);
    // 40 rad/s^2 at 100 mm is 4 m/s^2 — nearly half a g, from tipping alone.
    ASSERT_TRUE(minus_y < kG - 3.0);
}

// Coriolis: a ball ROLLING on a tilting plate is held by a different force
// than one sitting still on it, at the same instant and the same place.
void test_a_rolling_ball_is_held_differently_from_a_still_one() {
    const TableKinematics tk(plate());

    PlateMotion m = still(tk, 0, 0, 0.2121);
    m.omega = Eigen::Vector3d(1.5, 0.0, 0.0);   // tilting about x

    const Eigen::Vector3d s(0.05, 0.0, kR);
    const double at_rest = normalAccel(m, s, Eigen::Vector2d(0, 0), kG);
    const double rolling = normalAccel(m, s, Eigen::Vector2d(0, 0.4), kG);
    ASSERT_TRUE(std::abs(rolling - at_rest) > 0.1);
    // Rolling along +y under a +x rotation is carried downward, so it presses
    // less hard.  Reverse the roll and it presses harder by the same amount.
    const double other_way = normalAccel(m, s, Eigen::Vector2d(0, -0.4), kG);
    ASSERT_NEAR(rolling + other_way, 2.0 * at_rest, 1e-9);
}

// ---------------------------------------------------------------------------
// Leaving, flying, and landing
// ---------------------------------------------------------------------------

RollingBallDynamics roller() {
    const TableParams tp = plate();
    return RollingBallDynamics(kPlateBall, PlateParams{tp.R_table, kG});
}

// A ball on a plate that is behaving itself never leaves, and the rolling path
// is bit-for-bit the one `stepBall` gives — the contact layer adds a question,
// not a different answer.
//
// `stepBall` has to be told the same normal force the contact layer computed,
// which since #23 is the plate's cosine share of gravity and not `g`.  On a
// motionless plate `normalAccel` reduces to exactly `quasiStaticNormalAccel`,
// so this stays an equality rather than becoming a tolerance.  Handing the bare
// roller `g` here instead would compare the contact model against a ball on a
// LEVEL plate, which at 2 deg of tilt disagrees in the fourth decimal — the
// divergence this assertion exists to catch.
void test_a_settled_plate_never_lets_go() {
    const TableKinematics tk(plate());
    const RollingBallDynamics dyn = roller();
    const TablePose pose{2.0 * kDeg, -1.0 * kDeg, 0.2121};
    const PlateMotion m = still(tk, pose.phi, pose.theta, pose.z_c);

    BallState b;
    b.rolling << 0.05, -0.02, 0.0, 0.0;
    Eigen::Vector4d bare = b.rolling;

    for (int k = 0; k < 240; ++k) {
        b = stepBallContact(dyn, b, m, m, pose, kR, kG, 1.0 / 60.0);
        bare = stepBall(dyn, bare, pose, quasiStaticNormalAccel(m, kG), 1.0 / 60.0);
        ASSERT_TRUE(!b.airborne);
    }
    for (int i = 0; i < 4; ++i) ASSERT_NEAR(b.rolling(i), bare(i), 1e-15);
}

// Drop the plate out from under it and the ball goes ballistic — and the arc
// it follows is the one gravity alone gives, checked against the closed form
// rather than against itself.
//
// The plate has to actually recede, not merely be given a downward velocity:
// a ball in free fall drops at g, so a plate that wants to leave it behind has
// to drop faster.  At 3g it does, and the gap opens.
void test_a_dropped_plate_launches_the_ball_on_a_parabola() {
    const TableKinematics tk(plate());
    const RollingBallDynamics dyn = roller();
    const double dt = 1.0 / 60.0;
    const double a_plate = -3.0 * kG;

    double zc = 0.2121, vzc = 0.0;
    PlateMotion prev = still(tk, 0, 0, zc);

    BallState b;
    b.rolling << 0.0, 0.0, 0.30, 0.0;              // rolling along +x

    Eigen::Vector3d p0, v0;
    int launched_at = -1;
    for (int k = 0; k < 40; ++k) {
        vzc += a_plate * dt;
        zc += vzc * dt;
        const TablePose pose{0.0, 0.0, zc};
        PlateMotion now = still(tk, 0, 0, zc);
        now.c_dot = Eigen::Vector3d(0, 0, vzc);
        now.c_ddot = Eigen::Vector3d(0, 0, a_plate);

        const bool was_rolling = !b.airborne;
        b = stepBallContact(dyn, b, now, prev, pose, kR, kG, dt);
        prev = now;

        if (was_rolling && b.airborne) {
            launched_at = k;
            // It left carrying its rolling speed, and the plate's motion.
            ASSERT_NEAR(b.flight_v(0), 0.30, 1e-12);
            ASSERT_TRUE(b.flight_v(2) < 0.0);
            // Rewind the one step it took after leaving, to get the launch.
            v0 = b.flight_v - Eigen::Vector3d(0, 0, -kG) * dt;
            p0 = b.flight_p - v0 * dt - 0.5 * Eigen::Vector3d(0, 0, -kG) * dt * dt;
        }
    }
    ASSERT_TRUE(launched_at == 0);      // the very first frame, at 3g
    ASSERT_TRUE(b.airborne);            // and the plate never caught it again

    // s = ut + at^2/2, in the world, from the launch.
    const double t = (40 - launched_at) * dt;
    ASSERT_NEAR(b.flight_p(0), p0(0) + v0(0) * t, 1e-12);
    ASSERT_NEAR(b.flight_p(2), p0(2) + v0(2) * t - 0.5 * kG * t * t, 1e-12);
}

// And it comes back.  Landing is inelastic: the approach speed along the
// normal is absorbed, the sideways carry is kept, and the ball is rolling
// again from where it touched down rather than from where it took off.
void test_the_ball_lands_and_rolls_on() {
    const TableKinematics tk(plate());
    const RollingBallDynamics dyn = roller();
    const TablePose pose{0.0, 0.0, 0.2121};
    const double dt = 1.0 / 60.0;
    const PlateMotion m = still(tk, 0, 0, pose.z_c);

    BallState b;
    b.airborne = true;
    b.flight_p = Eigen::Vector3d(0.04, 0.0, pose.z_c + kR + 0.05);  // 50 mm up
    b.flight_v = Eigen::Vector3d(0.20, 0.0, 0.0);                   // drifting +x

    int frames = 0;
    while (b.airborne && frames < 600) {
        b = stepBallContact(dyn, b, m, m, pose, kR, kG, dt);
        ++frames;
    }
    ASSERT_TRUE(!b.airborne);
    // Free fall from 50 mm is about 0.10 s; one frame either side is fine.
    ASSERT_NEAR(frames * dt, std::sqrt(2.0 * 0.05 / kG), 0.02);
    // Carried downrange while it fell, and still moving that way.
    ASSERT_TRUE(b.rolling(0) > 0.04);
    ASSERT_NEAR(b.rolling(2), 0.20, 1e-9);
    // The drop is gone, not turned into a bounce.
    ASSERT_NEAR(b.rolling(3), 0.0, 1e-12);
}

// The plate-frame view is the same physical ball in both phases.  Converting a
// rolling ball out to the world and straight back must return exactly what it
// started as — including the velocity, which is the part that can go wrong:
// the plate frame is moving, so the world velocity of a ball rolling on a
// tilting plate is not its rolling velocity, and the way back has to remove
// again precisely what the way out added.
void test_the_two_frames_describe_the_same_ball() {
    const TableKinematics tk(plate());

    // A plate that is doing everything at once: tilted, heaving, and rotating.
    PlateMotion m = still(tk, 4.0 * kDeg, -3.0 * kDeg, 0.2121);
    m.c_dot = Eigen::Vector3d(0.0, 0.0, 0.35);
    m.omega = Eigen::Vector3d(0.8, -0.5, 0.2);

    BallState rolling;
    rolling.rolling << 0.06, -0.03, 0.10, 0.05;

    BallState flying;
    flying.airborne = true;
    worldOf(rolling, m, kR, &flying.flight_p, &flying.flight_v);

    const Eigen::Matrix<double, 6, 1> a = plateFrame(rolling, m, kR);
    const Eigen::Matrix<double, 6, 1> b = plateFrame(flying, m, kR);
    for (int i = 0; i < 6; ++i) ASSERT_NEAR(b(i), a(i), 1e-12);

    // And the world velocity really is more than the rolling velocity — if it
    // were not, the round trip above would be proving nothing.
    ASSERT_TRUE((flying.flight_v - m.R * Eigen::Vector3d(0.10, 0.05, 0.0)).norm() > 0.1);
}

}  // namespace

// ---------------------------------------------------------------------------
// The trust gate
// ---------------------------------------------------------------------------
//
// A separation is a claim about the plate's VELOCITY, and those rates come
// through `-J_pose^-1 J_alpha`, which near a singularity amplifies without
// bound.  `rates_trustworthy` is what stops the contact model acting on rates
// the mechanism's own condition number says are meaningless — the guard
// CONTEXT.md gives a pull-quote to, and which nothing tested.

// A plate diving hard enough to leave the ball behind — but whose rates nobody
// believes — keeps the ball.  Refusing to act is the whole point of the guard:
// the alternative is a hop that is arithmetic rather than physics.
void test_untrustworthy_rates_never_separate_the_ball() {
    const TableKinematics tk(plate());
    const RollingBallDynamics dyn = roller();
    const double dt = 1.0 / 60.0;

    // Falling at 3g, which for a trusted plate is a launch (see the parabola
    // test above), and the only difference here is the flag.
    PlateMotion prev = still(tk, 0, 0, 0.2121);
    PlateMotion now = still(tk, 0, 0, 0.2121 - 0.5 * 3.0 * kG * dt * dt);
    now.c_dot = Eigen::Vector3d(0, 0, -3.0 * kG * dt);
    now.c_ddot = Eigen::Vector3d(0, 0, -3.0 * kG);

    BallState trusted;
    trusted.rolling << 0.03, 0.0, 0.0, 0.0;
    BallState doubted = trusted;

    now.rates_trustworthy = true;
    prev.rates_trustworthy = true;
    trusted = stepBallContact(dyn, trusted, now, prev, TablePose{0, 0, 0.2121},
                              kR, kG, dt);
    ASSERT_TRUE(trusted.airborne);

    now.rates_trustworthy = false;
    doubted = stepBallContact(dyn, doubted, now, prev, TablePose{0, 0, 0.2121},
                              kR, kG, dt);
    ASSERT_TRUE(!doubted.airborne);
}

// And the half of the gate that is easy to drop.  The accelerations are
// analytic in the servo lag now, but not entirely differenceless: `plateMotion`
// still differences `J_v` and `A` across the two frames, and `J_v` is
// `-J_pose^-1 J_alpha` — the very quantity that blows up near a singularity.
// So the FIRST believable frame after a singular stretch still carries a term
// built from an unbelievable one, and trust has to cover both ends or the guard
// leaks exactly the hop it exists to refuse.
void test_a_trustworthy_frame_after_an_untrusted_one_still_declines() {
    const TableKinematics tk(plate());
    const RollingBallDynamics dyn = roller();
    const double dt = 1.0 / 60.0;

    // A plate whose own acceleration says the ball should leave, reached on the
    // first frame out of a stretch nobody believed.
    PlateMotion prev = still(tk, 0, 0, 0.2121);
    prev.rates_trustworthy = false;

    PlateMotion now = still(tk, 0, 0, 0.2121);
    now.c_ddot = Eigen::Vector3d(0, 0, -3.0 * kG);
    now.rates_trustworthy = true;

    BallState b;
    b.rolling << 0.03, 0.0, 0.0, 0.0;

    // The plate alone would say separation; the gate says no.
    const double N = normalAccel(now, Eigen::Vector3d(0.03, 0.0, kR),
                                 Eigen::Vector2d::Zero(), kG);
    ASSERT_TRUE(N < 0.0);

    b = stepBallContact(dyn, b, now, prev, TablePose{0, 0, 0.2121}, kR, kG, dt);
    ASSERT_TRUE(!b.airborne);
}

// ---------------------------------------------------------------------------
// The accelerations, and why they are not a difference
// ---------------------------------------------------------------------------

// The bug this replaced, stated as the thing that must not happen again.
//
// `omega_dot` and `c_ddot` used to come from differencing `omega` and `c_dot`
// across two frames.  Those rates are `(cmd - alpha) / tau` carried through the
// Jacobian, and `cmd` is a zero-order hold: it STEPS whenever the loop changes
// its mind.  Differencing a step gives `1/dt`, so the estimator answered a
// command jump with a delta function — 101 rad/s^2 at the corner of a square
// path on a thirty-second lap, where the setpoint is crawling at 24 mm/s.
//
// The analytic answer for the same jump is `alpha_ddot = -alpha_dot / tau`,
// which is what a servo with a 0.05 s lag actually does, and it does not depend
// on `dt` at all.  So: step the command, and check that the acceleration the
// plate reports is the servo's and not the frame rate's.
void test_a_stepped_command_is_not_an_infinite_acceleration() {
    const TableKinematics tk(plate());
    const double tau = 0.05;

    // The legs at home, and a command 5 degrees away on one of them — the size
    // of jump a corner produces.  Nothing has moved yet, so this is the frame
    // the step arrives on.
    const double home = M_PI / 4.0;
    const std::array<double, 3> alpha = {home, home, home};
    const std::array<double, 3> cmd = {home + 5.0 * kDeg, home, home};
    std::array<double, 3> adot{};
    for (int i = 0; i < 3; ++i) adot[i] = (cmd[i] - alpha[i]) / tau;

    const TablePose pose = tk.home_pose(home);
    const std::array<double, 3> addot = servoAccel(adot, tau);
    for (int i = 0; i < 3; ++i) ASSERT_NEAR(addot[i], -adot[i] / tau, 1e-15);

    // Two frame rates, and the acceleration is the same at both.  A difference
    // could not manage that: it would double when `dt` halved.
    const PlateMotion prev = plateMotion(tk, pose, alpha, {0.0, 0.0, 0.0});
    const PlateMotion slow =
        plateMotion(tk, pose, alpha, adot, addot, &prev, 1.0 / 60.0);
    const PlateMotion fast =
        plateMotion(tk, pose, alpha, adot, addot, &prev, 1.0 / 240.0);
    ASSERT_NEAR(slow.omega_dot.norm(), fast.omega_dot.norm(), 1e-9);
    ASSERT_NEAR(slow.c_ddot(2), fast.c_ddot(2), 1e-9);

    // And it is the servo's own number.  With the legs still at home, the
    // acceleration is `J_v alpha_ddot = -J_v alpha_dot / tau`, which is exactly
    // `-omega / tau`: 0.82 rad/s of plate rate through a 0.05 s lag is
    // 16.5 rad/s^2, and the table drops at 2.47 m/s^2 while it tilts.
    ASSERT_NEAR(slow.omega_dot.norm(), slow.omega.norm() / tau, 1e-9);
    ASSERT_NEAR(slow.omega_dot.norm(), 16.455, 1e-3);
    ASSERT_NEAR(slow.c_ddot(2), -2.468, 1e-3);

    // **And the size of the old error, exactly.**  On the frame a command
    // steps, `omega` goes from nothing to its full value, so differencing it
    // reports `|omega| / dt` against the true `|omega| / tau` — an overstatement
    // by `tau / dt`, which is 3x at 60 Hz and 12x at 240.  It is not noise that
    // averages out and it does not shrink with a finer step; it grows.
    ASSERT_NEAR(slow.omega.norm() / (1.0 / 60.0),
                slow.omega_dot.norm() * tau * 60.0, 1e-9);
    ASSERT_NEAR(slow.omega.norm() / (1.0 / 240.0),
                slow.omega_dot.norm() * tau * 240.0, 1e-9);
}

// The analytic answer, checked against a finite difference where a finite
// difference is actually trustworthy.
//
// The whole argument for going analytic is that the leg command steps, and a
// difference cannot differentiate a step.  Drive the legs with a SMOOTH
// sinusoid instead — no steps anywhere — and the two must agree, or the
// analytic expression is simply wrong.
//
// It also shows the two differenced terms are not optional.  `J_v` and `A` are
// themselves changing as the plate moves, and dropping their rates gets
// `omega_dot` wrong by 13% and `c_ddot` wrong in SIGN.
void test_the_analytic_acceleration_matches_a_difference_where_one_is_valid() {
    const TableKinematics tk(plate());
    const double home = M_PI / 4.0, w = 3.0, amp = 6.0 * kDeg;
    auto leg      = [&](double t, int i) { return home + amp * std::sin(w * t + i * 2.094); };
    auto leg_dot  = [&](double t, int i) { return amp * w * std::cos(w * t + i * 2.094); };
    auto leg_ddot = [&](double t, int i) { return -amp * w * w * std::sin(w * t + i * 2.094); };

    // March the pose seed forward so `solve_pose` stays on the right branch —
    // the constraint equations have a second root with the table folded flat
    // (#22), and a cold seed at t0 can land on it.
    TablePose seed = tk.home_pose(home);
    auto motion_at = [&](double t, const PlateMotion* prev, double dt) {
        const std::array<double, 3> a{leg(t, 0), leg(t, 1), leg(t, 2)};
        const std::array<double, 3> ad{leg_dot(t, 0), leg_dot(t, 1), leg_dot(t, 2)};
        const std::array<double, 3> add{leg_ddot(t, 0), leg_ddot(t, 1), leg_ddot(t, 2)};
        const FKResult fk = tk.solve_pose(a, seed);
        seed = fk.pose;
        return plateMotion(tk, fk.pose, a, ad, add, prev, dt);
    };
    for (double t = 0.0; t < 0.5; t += 1.0 / 600.0) motion_at(t, nullptr, 0.0);

    const double t0 = 0.5, h = 1e-5;
    const TablePose at_t0 = seed;
    const PlateMotion before = (seed = at_t0, motion_at(t0 - h, nullptr, 0.0));
    const PlateMotion after  = (seed = at_t0, motion_at(t0 + h, nullptr, 0.0));
    const Eigen::Vector3d omega_dot_fd = (after.omega - before.omega) / (2.0 * h);
    const Eigen::Vector3d c_ddot_fd = (after.c_dot - before.c_dot) / (2.0 * h);

    seed = at_t0;
    const PlateMotion prev = motion_at(t0 - 1.0 / 600.0, nullptr, 0.0);
    seed = at_t0;
    const PlateMotion now = motion_at(t0, &prev, 1.0 / 600.0);

    for (int i = 0; i < 3; ++i) ASSERT_NEAR(now.omega_dot(i), omega_dot_fd(i), 2e-3);
    ASSERT_NEAR(now.c_ddot(2), c_ddot_fd(2), 2e-3);

    // And the same numbers with `prev` withheld, which drops the `J_v` and `A`
    // rate terms: 0.617 against the true 0.546, and a `c_ddot` of +0.017 where
    // the truth is -0.006.  Small absolute numbers, but the sign of the heave is
    // the sign of whether the plate is pressing the ball or leaving it.
    seed = at_t0;
    const PlateMotion partial = motion_at(t0, nullptr, 0.0);
    ASSERT_TRUE(std::abs(partial.omega_dot(1) - omega_dot_fd(1)) > 0.05);
    ASSERT_TRUE(partial.c_ddot(2) * c_ddot_fd(2) < 0.0);
}

// A plate whose legs are moving at a CONSTANT rate has no acceleration to
// speak of — `alpha_ddot = -alpha_dot / tau` is what a lag does when it is
// chasing, and a leg already at its commanded angle is not.
void test_legs_at_their_command_produce_no_acceleration() {
    const TableKinematics tk(plate());
    const double home = M_PI / 4.0;
    const std::array<double, 3> alpha = {home, home, home};
    const TablePose pose = tk.home_pose(home);

    const PlateMotion m = plateMotion(tk, pose, alpha, {0.0, 0.0, 0.0},
                                      servoAccel({0.0, 0.0, 0.0}, 0.05));
    ASSERT_NEAR(m.omega.norm(), 0.0, 1e-12);
    ASSERT_NEAR(m.omega_dot.norm(), 0.0, 1e-12);
    ASSERT_NEAR(m.c_ddot.norm(), 0.0, 1e-12);
    ASSERT_NEAR(normalAccel(m, Eigen::Vector3d(0.05, 0.0, kR),
                            Eigen::Vector2d::Zero(), 9.81),
                9.81, 1e-9);
}

// The threshold is the application's own "Poor" line, not a second opinion.
void test_the_trust_threshold_is_the_condition_number_the_app_shows() {
    ASSERT_NEAR(kRatesUntrustworthyAbove, 20.0, 1e-15);
}

int main() {
    test_a_still_level_plate_holds_the_ball_at_one_g();
    test_a_tilted_still_plate_holds_the_cosine_share();
    test_heave_cancels_the_normal_force_at_one_g();
    test_steady_rotation_pulls_on_the_balls_own_radius();
    test_angular_acceleration_lifts_one_side_and_drops_the_other();
    test_a_rolling_ball_is_held_differently_from_a_still_one();
    test_a_settled_plate_never_lets_go();
    test_a_dropped_plate_launches_the_ball_on_a_parabola();
    test_the_ball_lands_and_rolls_on();
    test_the_two_frames_describe_the_same_ball();
    test_untrustworthy_rates_never_separate_the_ball();
    test_a_trustworthy_frame_after_an_untrusted_one_still_declines();
    test_the_trust_threshold_is_the_condition_number_the_app_shows();
    test_a_stepped_command_is_not_an_infinite_acceleration();
    test_legs_at_their_command_produce_no_acceleration();
    test_the_analytic_acceleration_matches_a_difference_where_one_is_valid();
    std::printf("test_ball_contact: all passed\n");
    return 0;
}
