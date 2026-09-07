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

// And it comes back — and bounces.  The first arrival reflects the approach
// speed at `kRestitution`, the sideways carry is untouched, and the ball is
// still in the air afterwards rather than stuck to the plate.
void test_the_ball_bounces_when_it_arrives() {
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
    while (frames < 600) {
        const double before = b.flight_v(2);
        b = stepBallContact(dyn, b, m, m, pose, kR, kG, dt);
        ++frames;
        if (b.flight_v(2) > 0.0 && before < 0.0) break;
    }
    ASSERT_TRUE(b.airborne);
    // Free fall from 50 mm is about 0.10 s; one frame either side is fine.
    ASSERT_NEAR(frames * dt, std::sqrt(2.0 * 0.05 / kG), 0.02);
    // It left at e times the speed a 50 mm drop arrives at.  Within a frame's
    // worth of gravity, because the impact is resolved inside the frame by
    // interpolation rather than exactly — the exact arithmetic is pinned by
    // `test_a_rising_plate_throws_the_ball_harder`, which arranges an impact at
    // a known instant.
    ASSERT_NEAR(b.flight_v(2), kRestitution * std::sqrt(2.0 * kG * 0.05),
                kG * dt);
    // Carried downrange while it fell, and the sideways carry is untouched by
    // the impact: no tangential impulse, deliberately.
    ASSERT_TRUE(b.flight_p(0) > 0.04);
    ASSERT_NEAR(b.flight_v(0), 0.20, 1e-12);
    // And it ends the frame ABOVE the surface, having flown out the remainder
    // of the frame under the rebound, so the next frame starts a flight rather
    // than finding itself still touching.
    ASSERT_TRUE(plateFrame(b, m, kR)(2) > kR);
}

// The bounce train ends, and it ends by running out of rebound rather than by
// running out of patience.  When it does, the ball rolls on from where it
// touched down carrying the speed it had sideways.
void test_the_bounce_train_terminates() {
    const TableKinematics tk(plate());
    const RollingBallDynamics dyn = roller();
    const TablePose pose{0.0, 0.0, 0.2121};
    const double dt = 1.0 / 60.0;
    const PlateMotion m = still(tk, 0, 0, pose.z_c);

    BallState b;
    b.airborne = true;
    b.flight_p = Eigen::Vector3d(0.04, 0.0, pose.z_c + kR + 0.05);
    b.flight_v = Eigen::Vector3d(0.20, 0.0, 0.0);

    int frames = 0, bounces = 0;
    while (b.airborne && frames < 6000) {
        const double before = b.flight_v(2);
        b = stepBallContact(dyn, b, m, m, pose, kR, kG, dt);
        ++frames;
        if (b.airborne && b.flight_v(2) > 0.0 && before < 0.0) ++bounces;
        // No frame ever ends with the ball inside the plate.  This is what the
        // sub-frame impact buys: bouncing at the frame boundary would leave the
        // ball below the surface with the whole frame's extra fall counted as
        // approach speed, which is how the train came to gain energy instead of
        // losing it.
        if (b.airborne) ASSERT_TRUE(plateFrame(b, m, kR)(2) >= kR);
    }
    ASSERT_TRUE(!b.airborne);

    // e = 0.94 keeps 88% of the energy per bounce, so the train is long: from
    // 0.99 m/s it takes ln(g dt / 2 / u0) / ln(e) = 40 rebounds to fall under
    // the floor, and about 2 u0 e / (g (1 - e)) = 3.0 s of bouncing.
    ASSERT_TRUE(bounces > 30 && bounces < 50);
    ASSERT_TRUE(frames * dt > 2.0 && frames * dt < 4.0);

    // Still drifting the way it was, still on the plate.
    ASSERT_NEAR(b.rolling(2), 0.20, 1e-9);
    ASSERT_TRUE(b.rolling(0) > 0.04);
}

// **Restitution acts on the RELATIVE normal velocity at the contact point**,
// which is the whole difference between a bounce model and a bounce decoration:
// a plate driven up into a falling ball throws it harder than it arrived, and a
// plate running away catches it softly.  The same term is what makes a
// deliberate hop possible at all.
void test_a_rising_plate_throws_the_ball_harder() {
    const TableKinematics tk(plate());
    const RollingBallDynamics dyn = roller();
    const TablePose pose{0.0, 0.0, 0.2121};
    const double dt = 1.0 / 60.0;
    const double drop = 1.0;   // [m/s] arriving

    // Three plates at the same place, differing only in how they are moving.
    // Flush with the surface and already moving into it, so the impact is at
    // the frame's open and the arrival speed is exactly `drop`.
    auto arrive = [&](double plate_vz) {
        PlateMotion m = still(tk, 0, 0, pose.z_c);
        m.c_dot = Eigen::Vector3d(0, 0, plate_vz);
        BallState b;
        b.airborne = true;
        b.flight_p = Eigen::Vector3d(0.0, 0.0, pose.z_c + kR);
        b.flight_v = Eigen::Vector3d(0.0, 0.0, -drop);
        b = stepBallContact(dyn, b, m, m, pose, kR, kG, dt);
        return b;
    };

    const BallState still_plate = arrive(0.0);
    const BallState rising = arrive(0.5);
    const BallState receding = arrive(-0.5);

    ASSERT_TRUE(still_plate.airborne && rising.airborne);

    // The approach is (v_ball - v_plate) . n, so a plate rising at 0.5 m/s is
    // approached 0.5 m/s faster and the rebound is e times that much more —
    // and the ball leaves in the WORLD carrying the plate's velocity too.  Each
    // then flies out the remainder of the frame, which is all of it.
    auto leaves_at = [&](double plate_vz) {
        return plate_vz + kRestitution * (drop + plate_vz) - kG * dt;
    };
    ASSERT_NEAR(still_plate.flight_v(2), leaves_at(0.0), 1e-12);
    ASSERT_NEAR(rising.flight_v(2), leaves_at(0.5), 1e-12);
    ASSERT_TRUE(rising.flight_v(2) > still_plate.flight_v(2) + 0.9);

    // And a plate running away at 0.5 m/s takes half the arrival out of the
    // rebound, which is the same arithmetic pointed the other way.
    ASSERT_NEAR(receding.flight_v(2), leaves_at(-0.5), 1e-12);
    ASSERT_TRUE(receding.flight_v(2) < still_plate.flight_v(2) - 0.9);
}

// The floor is derived from the frame rate, not chosen against the passive
// case: a rebound whose whole flight fits inside one frame is not small, it is
// unrepresentable.  So it moves when `dt` does, and a ball that is genuinely
// hopping cannot reach it — the shipped tuning's own passive hop leaves the
// plate some twenty times faster.
void test_the_bounce_floor_is_the_frame_rate_and_nothing_else() {
    ASSERT_NEAR(bounceFloorSpeed(kG, 1.0 / 60.0), 0.5 * kG / 60.0, 1e-15);
    ASSERT_NEAR(bounceFloorSpeed(kG, 1.0 / 60.0), 0.08175, 1e-5);
    // Halve the step and the simulation resolves twice as fine a bounce.
    ASSERT_NEAR(bounceFloorSpeed(kG, 1.0 / 120.0),
                0.5 * bounceFloorSpeed(kG, 1.0 / 60.0), 1e-15);

    // A rebound exactly at the floor rises for exactly one frame: u dt - g dt^2
    // / 2 = 0 at u = g dt / 2.  That is the sense in which it is the smallest
    // representable bounce and not a threshold anybody picked.
    const double dt = 1.0 / 60.0;
    const double u = bounceFloorSpeed(kG, dt);
    ASSERT_NEAR(u * dt - 0.5 * kG * dt * dt, 0.0, 1e-18);
}

// Either side of the floor, through `stepBallContact` rather than in
// arithmetic: just above it the ball bounces, just below it the ball is landed
// and rolling.
void test_a_rebound_below_the_floor_lands_instead() {
    const TableKinematics tk(plate());
    const RollingBallDynamics dyn = roller();
    const TablePose pose{0.0, 0.0, 0.2121};
    const double dt = 1.0 / 60.0;
    const PlateMotion m = still(tk, 0, 0, pose.z_c);
    const double floor_u = bounceFloorSpeed(kG, dt);

    // Arrive at a chosen speed, by starting flush with the surface so the one
    // step both applies gravity and makes contact.
    auto arrive_at = [&](double approach) {
        BallState b;
        b.airborne = true;
        b.flight_p = Eigen::Vector3d(0.0, 0.0, pose.z_c + kR);
        b.flight_v = Eigen::Vector3d(0.0, 0.0, -approach);
        return stepBallContact(dyn, b, m, m, pose, kR, kG, dt);
    };

    ASSERT_TRUE(arrive_at(1.02 * floor_u / kRestitution).airborne);
    ASSERT_TRUE(!arrive_at(0.98 * floor_u / kRestitution).airborne);
}

// A bounce is a claim about the plate's velocity, so it is refused on the same
// grounds a separation is.  Near a singularity the rates are arithmetic rather
// than physics, and a rebound built on them would be a launch the mechanism
// never performed — the ball lands instead, which is what it did before this
// model could bounce at all.
void test_untrusted_rates_land_the_ball_rather_than_bouncing_it() {
    const TableKinematics tk(plate());
    const RollingBallDynamics dyn = roller();
    const TablePose pose{0.0, 0.0, 0.2121};
    const double dt = 1.0 / 60.0;

    PlateMotion m = still(tk, 0, 0, pose.z_c);
    m.c_dot = Eigen::Vector3d(0, 0, 12.0);   // the wild rate a bad solve returns

    BallState b;
    b.airborne = true;
    b.flight_p = Eigen::Vector3d(0.0, 0.0, pose.z_c + kR);
    b.flight_v = Eigen::Vector3d(0.0, 0.0, -0.5);

    m.rates_trustworthy = true;
    const BallState trusted = stepBallContact(dyn, b, m, m, pose, kR, kG, dt);
    ASSERT_TRUE(trusted.airborne);
    ASSERT_TRUE(trusted.flight_v(2) > 10.0);   // launched, by the arithmetic

    m.rates_trustworthy = false;
    const BallState doubted = stepBallContact(dyn, b, m, m, pose, kR, kG, dt);
    ASSERT_TRUE(!doubted.airborne);
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
    test_the_ball_bounces_when_it_arrives();
    test_the_bounce_train_terminates();
    test_a_rising_plate_throws_the_ball_harder();
    test_the_bounce_floor_is_the_frame_rate_and_nothing_else();
    test_a_rebound_below_the_floor_lands_instead();
    test_untrusted_rates_land_the_ball_rather_than_bouncing_it();
    test_the_two_frames_describe_the_same_ball();
    test_untrustworthy_rates_never_separate_the_ball();
    test_a_trustworthy_frame_after_an_untrusted_one_still_declines();
    test_the_trust_threshold_is_the_condition_number_the_app_shows();
    std::printf("test_ball_contact: all passed\n");
    return 0;
}
