// tests/test_hop_drive.cpp
//
// The hopping controller: the plate throwing the ball on purpose while the
// loop goes on tracking a path.  Two layers, as everything about the contact
// is: the cycle on its own, driven by hand where the answer is known, and then
// the loop closed around the nonlinear plate through `stepSim`.
//
// **What is being tested is that the hop is #23's constraint with the sign
// released and nothing else.**  `holdContactDown` drives the contact point's
// normal rate to zero because `u_rebound = e u + p (1 + e)` makes `p <= 0` a
// proof; a throw is the same `heaveToContactRate` solve asked for `p > 0` at a
// contact the cycle chose.  So the interesting assertions are that one
// function serves both signs, that the flight is still the constraint's, and
// that the two gates on arming keep the ball on the plate.
//
// See issue #34.
#include "hop_drive.h"

#include "auto_balance.h"
#include "ball_contact.h"
#include "cascade_fixture.h"
#include "setpoint_path.h"
#include "sim_step.h"
#include "table_kinematics.h"
#include "test_helpers.h"

#include <cmath>
#include <cstdio>

using namespace caliburn;

namespace {

constexpr double kBallRadius = kFixtureBallRadius;
const double kHome = M_PI / 4.0;

// ---------------------------------------------------------------------------
// Layer 1: the cycle, with no plate under it
// ---------------------------------------------------------------------------
//
// `stepHop` knows one scalar about the mechanism — the plate's heave — and
// answers with one scalar about the contact point.  That is the whole of its
// interface, and it is why these tests need no kinematics.

/// A sizing whose numbers are round, for the cases where the arithmetic is the
/// thing under test rather than the plant.
HopSizing roundSizing() {
    HopSizing s;
    s.rise_rate = 1.0;      // m/s
    s.stroke_s = 0.1;       // s
    s.drop_m = 0.05;        // m, the area under the ramp
    s.flight_s = 2.0;
    s.max_ball_speed = 0.1;
    return s;
}

/// One frame in which everything the gates look at is in order.
HopFrame ready(double plate_z) {
    HopFrame f;
    f.enabled = true;
    f.loop_has_room = true;
    f.plate_z = plate_z;
    return f;
}

/// Drive the cycle with the plate doing exactly what it is asked, so the
/// travel conditions behave.  Returns the demand and moves `z`.
HopDemand run(const HopSizing& size, HopCycle& c, double& z, double dt) {
    const HopDemand d = stepHop(size, c, ready(z), dt);
    if (d.drive) z += d.contact_rate * dt;
    return d;
}

void test_the_hop_is_sized_by_the_servo_the_ball_and_the_rim() {
    const double g = 9.81, tau = 0.05, room = 0.28;
    const HopSizing s = hopSizing(g, tau, room);

    // The throw's rate is a margin over `g tau`, which is the rate a plate has
    // to be rising at before its own lag can decelerate it past gravity and
    // let the ball go.  Below that there is no hop to be had at any stroke.
    ASSERT_NEAR(s.rise_rate, kHopRiseMargin * g * tau, 1e-12);
    ASSERT_TRUE(s.rise_rate > g * tau);

    ASSERT_NEAR(s.stroke_s, kHopStrokeLags * tau, 1e-12);

    // The strokes ramp, so the travel is the area under a triangle.  Half of
    // `rate * time`, not all of it — see `HopSizing::drop_m`.
    ASSERT_NEAR(s.drop_m, 0.5 * s.rise_rate * s.stroke_s, 1e-12);

    // The first flight is `2 u / g` and the train after it is `1 / (1 - e)`
    // times that, from the one `kRestitution` the bounce itself uses.
    ASSERT_NEAR(s.flight_s, 2.0 * s.rise_rate / (g * (1.0 - kRestitution)), 1e-12);
    ASSERT_NEAR(s.max_ball_speed, room / s.flight_s, 1e-12);

    // The numbers this ticket's measurements were taken against.
    ASSERT_NEAR(s.rise_rate, 0.7358, 1e-3);
    ASSERT_NEAR(s.drop_m, 0.0368, 1e-4);
    ASSERT_NEAR(s.flight_s, 2.50, 0.01);
    ASSERT_NEAR(s.max_ball_speed, 0.112, 1e-3);

    // A smaller plate refuses to hop a ball that is moving as fast, because a
    // throw buys the same flight and there is less room to come down in.
    ASSERT_TRUE(hopSizing(g, tau, 0.5 * room).max_ball_speed <
                0.51 * s.max_ball_speed);
}

// Switched off, the cycle is not paused — it is forgotten.  Resuming into a
// half-charged plate would leave the legs 37 mm low with nothing to spend it on.
void test_disabling_the_hop_forgets_the_cycle() {
    const HopSizing size = roundSizing();
    HopCycle c;
    double z = 0.2;
    for (int k = 0; k < 4; ++k) run(size, c, z, 1.0 / 60.0);
    ASSERT_TRUE(c.phase == HopPhase::Charge);

    const HopDemand off = stepHop(size, c, HopFrame{}, 1.0 / 60.0);
    ASSERT_TRUE(!off.drive);
    ASSERT_TRUE(c.phase == HopPhase::Off);
    ASSERT_NEAR(c.phase_s, 0.0, 1e-12);
}

// Down together, then up: the shape the ticket asks for, and the shape that
// makes the throw a throw rather than a lift.
void test_a_hop_charges_down_before_it_throws_up() {
    const HopSizing size = roundSizing();
    HopCycle c;
    double z = 0.2;
    const double dt = 1.0 / 60.0;

    bool saw_charge = false, saw_throw = false;
    double lowest = z;
    for (int k = 0; k < 40; ++k) {
        const HopDemand d = run(size, c, z, dt);
        ASSERT_TRUE(d.drive);
        if (c.phase == HopPhase::Charge) {
            ASSERT_TRUE(d.contact_rate <= 0.0);
            ASSERT_TRUE(!saw_throw);      // charge first, and only once
            saw_charge = true;
        } else {
            ASSERT_TRUE(c.phase == HopPhase::Throw);
            ASSERT_TRUE(d.contact_rate >= 0.0);
            saw_throw = true;
        }
        lowest = std::min(lowest, z);
        if (saw_throw && z >= 0.2) break;
    }
    ASSERT_TRUE(saw_charge && saw_throw);

    // And the charge bought exactly what it said it would.
    ASSERT_NEAR(0.2 - lowest, size.drop_m, 2.0 * size.rise_rate * dt);
}

// **The ramp is the decision.**  A rate asked for flat arrives inside one
// frame, the plate steps past the release threshold, and what the ball gets is
// an impact rather than a carry — measured, a 53 mm hop instead of a 24 mm one,
// with the no-pumping constraint given up by 0.44 m/s.  So the demand has to
// grow across the stroke rather than appear at full size.
void test_the_strokes_ramp_rather_than_step() {
    const HopSizing size = roundSizing();
    HopCycle c;
    double z = 0.2;
    const double dt = size.stroke_s / 10.0;

    double last = 0.0;
    for (int k = 0; k < 5; ++k) {
        const HopDemand d = run(size, c, z, dt);
        ASSERT_TRUE(c.phase == HopPhase::Charge);
        const double asked = -d.contact_rate;
        ASSERT_TRUE(asked > last);                 // strictly growing
        ASSERT_TRUE(asked <= size.rise_rate);      // and never past the target
        last = asked;
    }
    // A tenth of the way in it is asking for a tenth of the rate, not all of it.
    HopCycle fresh;
    double z2 = 0.2;
    ASSERT_NEAR(-run(size, fresh, z2, dt).contact_rate, 0.1 * size.rise_rate, 1e-12);
}

// The strokes end on TRAVEL, and the timer is only the backstop.  A plate that
// descends faster than it was asked hands over early; one that cannot descend
// at all hands over on the window rather than charging for ever.
void test_the_strokes_end_on_travel_with_the_timer_as_a_backstop() {
    const HopSizing size = roundSizing();
    const double dt = 1.0 / 60.0;

    // Descending far faster than asked: the charge is over as soon as the
    // plate is `drop_m` down, well inside the stroke's window.
    HopCycle fast;
    double z = 0.2;
    int k = 0;
    for (; k < 100 && fast.phase != HopPhase::Throw; ++k) {
        stepHop(size, fast, ready(z), dt);
        z -= 0.02;                                  // 20 mm a frame
    }
    ASSERT_TRUE(fast.phase == HopPhase::Throw);
    ASSERT_TRUE(k * dt < size.stroke_s);

    // A plate that does not move at all: the window is what ends the charge,
    // so the throw still happens and the cycle does not stall.
    HopCycle stuck;
    double held = 0.2;
    int frames = 0;
    for (; frames < 100 && stuck.phase != HopPhase::Throw; ++frames)
        stepHop(size, stuck, ready(held), dt);
    ASSERT_TRUE(stuck.phase == HopPhase::Throw);
    ASSERT_NEAR(frames * dt, size.stroke_s, 2.0 * dt);
}

// **What makes `drop_m` a bound.**  A throw that fails to let go recharges, and
// before the arming heave was remembered each cycle started from wherever the
// last one left the plate — measured, 101 mm below home on a path whose stated
// cost was 37, and four balls lost of eighteen settings.
void test_a_failed_throw_does_not_march_the_plate_down() {
    const HopSizing size = roundSizing();
    HopCycle c;
    const double dt = 1.0 / 60.0;
    const double z0 = 0.2;
    double z = z0;

    double lowest = z0;
    for (int k = 0; k < 600; ++k) {     // ten cycles' worth, never letting go
        run(size, c, z, dt);
        lowest = std::min(lowest, z);
    }
    ASSERT_TRUE(z0 - lowest < size.drop_m + 2.0 * size.rise_rate * dt);
}

// The hop does not drive a flight.  From the frame the ball leaves to the frame
// it lands the command is `holdContactDown`'s, which is what bounds the hop by
// the throw that started it rather than by a sweep.
void test_the_flight_belongs_to_the_constraint() {
    const HopSizing size = roundSizing();
    HopCycle c;
    double z = 0.2;
    const double dt = 1.0 / 60.0;
    for (int k = 0; k < 8; ++k) run(size, c, z, dt);

    HopFrame flying = ready(z);
    flying.airborne = true;
    const HopDemand d = stepHop(size, c, flying, dt);
    ASSERT_TRUE(!d.drive);
    ASSERT_TRUE(c.phase == HopPhase::Flight);
    ASSERT_EQ(c.thrown, 1);

    // And a whole train of it changes nothing.
    for (int k = 0; k < 200; ++k) {
        ASSERT_TRUE(!stepHop(size, c, flying, dt).drive);
        ASSERT_TRUE(c.phase == HopPhase::Flight);
    }
    ASSERT_EQ(c.thrown, 1);

    // Landing arms the next one — once per landing, never inside the train.
    run(size, c, z, dt);
    ASSERT_TRUE(c.phase == HopPhase::Charge);
}

// Only a throw that let the ball go is counted.  A window that expires with the
// ball still on the plate is the failure this number exists to make visible.
void test_only_a_throw_that_let_go_is_counted() {
    const HopSizing size = roundSizing();
    HopCycle c;
    double z = 0.2;
    const double dt = 1.0 / 60.0;
    for (int k = 0; k < 600; ++k) run(size, c, z, dt);
    ASSERT_TRUE(c.thrown == 0);
}

// The two gates, and they gate ARMING rather than a stroke already under way.
void test_a_throw_is_refused_while_the_ball_is_fast_or_the_loop_is_on_its_stops() {
    const HopSizing size = roundSizing();
    const double dt = 1.0 / 60.0;
    double z = 0.2;

    HopCycle fast;
    HopFrame rushing = ready(z);
    rushing.ball_speed = size.max_ball_speed * 1.01;
    ASSERT_TRUE(!stepHop(size, fast, rushing, dt).drive);
    ASSERT_TRUE(fast.phase == HopPhase::Waiting);

    HopCycle pinned;
    HopFrame on_the_stops = ready(z);
    on_the_stops.loop_has_room = false;
    ASSERT_TRUE(!stepHop(size, pinned, on_the_stops, dt).drive);
    ASSERT_TRUE(pinned.phase == HopPhase::Waiting);

    // Just under the speed gate it arms.
    HopCycle ok;
    HopFrame just_slow = ready(z);
    just_slow.ball_speed = size.max_ball_speed * 0.99;
    ASSERT_TRUE(stepHop(size, ok, just_slow, dt).drive);

    // And a hop already charging is carried through to its throw: a plate left
    // low is worse than a plate that finished what it started.
    HopCycle busy;
    double zb = 0.2;
    for (int k = 0; k < 3; ++k) run(size, busy, zb, dt);
    ASSERT_TRUE(busy.phase == HopPhase::Charge);
    HopFrame both_gates_shut = ready(zb);
    both_gates_shut.loop_has_room = false;
    both_gates_shut.ball_speed = 10.0;
    ASSERT_TRUE(stepHop(size, busy, both_gates_shut, dt).drive);
}

// ---------------------------------------------------------------------------
// Layer 1b: one heave actuator, two signs
// ---------------------------------------------------------------------------

// **The whole of #34's claim to be #23's second half.**  `holdContactDown`
// drives `p` to zero and the throw drives it positive, and it is the same
// solve: the tilt rates come out bit for bit unchanged at every target.
void test_the_heave_solve_hits_its_target_at_either_sign() {
    const auto models = getBuiltinModels();
    const TableParams tp = cascadeMechanism(cascadeModel(models).params);
    const TableKinematics tk(tp);
    // Off the symmetric pose on purpose, where a uniform (1, 1, 1) nudge would
    // tilt the plate as it lifted it and this must not.
    const std::array<double, 3> alpha = {kHome + 0.02, kHome - 0.03, kHome + 0.01};
    const FKResult fk = tk.solve_pose(alpha, tk.home_pose(kHome));
    ASSERT_TRUE(fk.converged);

    AutoBalanceDesign d;
    d.home_leg_rad = kHome;
    d.servo_tau = 0.05;
    d.mechanism = tp;

    const std::array<double, 3> cmd = {kHome + 0.14, kHome + 0.02, kHome + 0.07};
    const PlateMotion before = movingPlate(tk, fk.pose, alpha, cmd, d.servo_tau);
    const Eigen::Vector3d s(0.04, -0.06, kBallRadius);

    for (double target : {-0.40, -0.10, 0.0, 0.10, 0.40}) {
        const std::array<double, 3> aimed =
            heaveToContactRate(tk, d, before, alpha, cmd, s, target);
        const PlateMotion after = movingPlate(tk, fk.pose, alpha, aimed, d.servo_tau);
        ASSERT_NEAR(contactNormalRate(after, s), target, 1e-9);
        // Tilt authority untouched, and exactly so.
        for (int i = 0; i < 3; ++i)
            ASSERT_NEAR(after.omega(i), before.omega(i), 1e-9);
    }

    // And the constraint IS the solve asked for zero, bit for bit.
    ASSERT_TRUE(commandedContactRate(tk, d, before, alpha, cmd, s) > 0.0);
    const std::array<double, 3> held = holdContactDown(tk, d, before, alpha, cmd, s);
    const std::array<double, 3> aimed =
        heaveToContactRate(tk, d, before, alpha, cmd, s, 0.0);
    for (int i = 0; i < 3; ++i) ASSERT_EQ(held[i], aimed[i]);
}

// **"All three legs down together, then up" is a claim, so it is checked.**
// The correction is `J_v^-1 (0, 0, dz)` rather than a uniform `(1, 1, 1)` —
// that is what leaves the tilt alone away from the symmetric pose — so the
// three increments are not equal.  What they are is the SAME SIGN and the same
// order of magnitude, which is what the words mean and what the panel tells a
// visitor to look for.  A correction that sent one leg the other way would
// still produce pure heave on paper and would not be a common-mode stroke.
void test_a_throw_moves_all_three_legs_the_same_way() {
    const auto models = getBuiltinModels();
    const TableParams tp = cascadeMechanism(cascadeModel(models).params);
    const TableKinematics tk(tp);
    const std::array<double, 3> alpha = {kHome + 0.02, kHome - 0.03, kHome + 0.01};
    const FKResult fk = tk.solve_pose(alpha, tk.home_pose(kHome));
    ASSERT_TRUE(fk.converged);

    AutoBalanceDesign d;
    d.home_leg_rad = kHome;
    d.servo_tau = 0.05;
    d.mechanism = tp;

    // A tracking command with a tilt in it, from an off-centre pose, so the
    // heave and the tilt are genuinely mixed in the triple.
    const std::array<double, 3> cmd = {kHome + 0.06, kHome - 0.01, kHome + 0.03};
    const PlateMotion plate = movingPlate(tk, fk.pose, alpha, cmd, d.servo_tau);
    const Eigen::Vector3d s(0.04, -0.06, kBallRadius);
    const double rest = commandedContactRate(tk, d, plate, alpha, cmd, s);

    for (double target : {rest - 0.40, rest + 0.40}) {
        const std::array<double, 3> aimed =
            heaveToContactRate(tk, d, plate, alpha, cmd, s, target);
        const double want = (target > rest) ? 1.0 : -1.0;
        double smallest = 1e9, largest = 0.0;
        for (int i = 0; i < 3; ++i) {
            const double step = aimed[i] - cmd[i];
            ASSERT_TRUE(step * want > 0.0);        // all three, the same way
            smallest = std::min(smallest, std::abs(step));
            largest = std::max(largest, std::abs(step));
        }
        // ...and together rather than one leg doing the work: the spread is
        // the plate's geometry, not a tilt smuggled in as heave.
        ASSERT_TRUE(largest < 2.0 * smallest);
    }
}

// The refusals are the constraint's, unchanged: no lag means no map from a
// command to a rate, and rates nobody vouches for are not rates to invert.
void test_the_heave_solve_declines_what_the_constraint_declines() {
    const auto models = getBuiltinModels();
    const TableParams tp = cascadeMechanism(cascadeModel(models).params);
    const TableKinematics tk(tp);
    const TablePose pose = tk.home_pose(kHome);
    const std::array<double, 3> alpha = {kHome, kHome, kHome};
    const Eigen::Vector3d s(0.05, 0.0, kBallRadius);

    AutoBalanceDesign d;
    d.home_leg_rad = kHome;
    d.servo_tau = 0.05;
    d.mechanism = tp;

    PlateMotion m = movingPlate(tk, pose, alpha, alpha, d.servo_tau);
    const std::array<double, 3> cmd = {kHome, kHome, kHome};

    m.rates_trustworthy = false;
    for (int i = 0; i < 3; ++i)
        ASSERT_EQ(heaveToContactRate(tk, d, m, alpha, cmd, s, 0.5)[i], cmd[i]);

    m.rates_trustworthy = true;
    AutoBalanceDesign lagless = d;
    lagless.servo_tau = 0.0;
    for (int i = 0; i < 3; ++i)
        ASSERT_EQ(heaveToContactRate(tk, lagless, m, alpha, cmd, s, 0.5)[i], cmd[i]);

    // With both in order it does move the command.
    ASSERT_TRUE(heaveToContactRate(tk, d, m, alpha, cmd, s, 0.5)[0] > cmd[0]);
}

// ---------------------------------------------------------------------------
// Layer 2: the loop, closed around the nonlinear plate
// ---------------------------------------------------------------------------

struct Hop {
    double mean_err = 0.0, max_err = 0.0;   ///< [m] after the first lap
    double max_radius = 0.0;                ///< [m] furthest the ball got out
    double apex = 0.0;                      ///< [m] highest above the surface
    double worst_rise = 0.0;                ///< [m/s] `p` in flight; #23's bound
    double lowest_z = 0.0;                  ///< [m] the heave the hop spent
    int airborne = 0, frames = 0, thrown = 0;
    int impacts_at_release = 0;             ///< throws that let go ON an impact
    double weakest_release = 1e9;           ///< [m/s] rebound at a release impact
    bool lost = false, changed_assembly = false;
};

/// The application's own step, with the hop switched on after one clean lap so
/// that what is measured is a hop rather than a start-up transient.
Hop follow(const ModelEntry& e, const Eigen::MatrixXd& K,
           const SetpointPath& asked, bool hop, double duration,
           double dt = 1.0 / 60.0) {
    const SimPlate plate = cascadePlate(e.params);
    SimInput in;
    in.dt = dt;
    in.design = cascadeDesign(e.params);
    in.design.K = K;
    in.path = plate.feasible(asked);

    const Eigen::Vector2d start = pathPoint(in.path, 0.0);
    SimState s = simStart(plate, in.design.home_leg_rad,
                          Eigen::Vector4d(start(0), start(1), 0.0, 0.0));

    Hop r;
    r.lowest_z = s.pose.z_c;
    double t = 0.0, sum = 0.0;
    int n = 0, thrown_before = 0;
    bool was_airborne = false;
    for (int k = 0; k < static_cast<int>(duration / dt); ++k) {
        in.hop_enabled = hop && t > in.path.period_s;
        const SimReport f = stepSim(plate, in, s);
        t += dt;
        ++r.frames;

        r.apex = std::max(r.apex, f.ball_plate(2) - SimPlate::ballRadius());
        r.max_radius = std::max(r.max_radius,
                                std::hypot(f.ball_plate(0), f.ball_plate(1)));
        r.lowest_z = std::min(r.lowest_z, s.pose.z_c);
        if (f.airborne) ++r.airborne;
        // `p` is only a constraint on a frame the ball ENTERED airborne — the
        // frame a separation happens on is one where the ball was on the plate
        // when the command was chosen.  See `SimReport::contact_normal_rate`.
        if (was_airborne && f.airborne)
            r.worst_rise = std::max(r.worst_rise, f.contact_normal_rate);
        // The frame a throw let go on.  A commanded hop is a SEPARATION, so the
        // bounce law — and with it the bounce floor — is not consulted there.
        if (!was_airborne && f.airborne && f.hops_thrown > thrown_before) {
            if (f.impact_approach != 0.0) {
                ++r.impacts_at_release;
                r.weakest_release = std::min(
                    r.weakest_release, -kRestitution * f.impact_approach);
            }
            thrown_before = f.hops_thrown;
        }
        was_airborne = f.airborne;

        r.thrown = f.hops_thrown;
        if (f.left_plate) r.lost = true;
        if (!onBuiltAssembly(plate.kinematics(), s.alpha_rad, s.pose))
            r.changed_assembly = true;
        if (t > in.path.period_s) {
            const double err = std::hypot(f.ball_plate(0) - f.setpoint(0),
                                          f.ball_plate(1) - f.setpoint(1));
            r.max_err = std::max(r.max_err, err);
            sum += err;
            ++n;
        }
    }
    r.mean_err = n ? sum / n : 0.0;
    return r;
}

SetpointPath openingCircle() {
    SetpointPath p;
    p.shape = PathShape::Circle;
    p.radius_m = 0.12;
    p.period_s = 10.0;
    return p;
}

// The headline: the plate throws the ball on purpose, repeatedly, while the
// ball goes on round the circle.
void test_the_plate_hops_the_ball_while_it_tracks() {
    const auto models = getBuiltinModels();
    const ModelEntry& e = cascadeModel(models);
    const Eigen::MatrixXd K = defaultGain(e);
    const SetpointPath p = openingCircle();

    const Hop flat = follow(e, K, p, false, 30.0);
    ASSERT_EQ(flat.airborne, 0);        // tracking alone never separates it
    ASSERT_EQ(flat.thrown, 0);

    const Hop hopping = follow(e, K, p, true, 30.0);

    // Repeatedly and deliberately: about one throw every two and a half
    // seconds of hopping, which is the train's own length.
    ASSERT_TRUE(hopping.thrown >= 6);
    ASSERT_TRUE(hopping.airborne > hopping.frames / 4);

    // **The hop height is bounded and stated.**  `rise_rate^2 / 2g` is 27.6 mm;
    // measured 24.7, because the ball is released a little under the rate the
    // ramp was heading for.
    ASSERT_NEAR(hopping.apex, 0.0247, 0.006);

    // And it still goes round the path.  3.6 mm of mean error without the hop,
    // 4.0 mm with it — the tilt is untouched by construction and what the hop
    // costs is the contact time the loop can steer in.
    ASSERT_NEAR(flat.mean_err, 0.0036, 0.0005);
    ASSERT_NEAR(hopping.mean_err, 0.0040, 0.0012);
    ASSERT_TRUE(hopping.max_err < 0.015);

    // The ball is never lost, and never near it: 130 mm out against a 280 mm
    // rim on a 120 mm path.
    ASSERT_TRUE(!hopping.lost);
    ASSERT_TRUE(!hopping.changed_assembly);
    ASSERT_TRUE(hopping.max_radius < 0.14);
}

// **The workspace cost, stated.**  A hop spends `drop_m` of heave and nothing
// else: the throw is cut short the frame the ball leaves and the flight is the
// constraint's.  What that costs is the tilt the plate can still be HELD at
// from a lower heave, which is `max_conditioned_tilt` — the same quantity the
// corner fillet is sized by.
void test_the_hop_spends_a_stated_amount_of_tilt_authority() {
    const auto models = getBuiltinModels();
    const ModelEntry& e = cascadeModel(models);
    const SimPlate plate = cascadePlate(e.params);
    const TableKinematics& tk = plate.kinematics();
    const double g = plate.gravity();
    const double home_z = tk.home_pose(cascadeHomeLegAngle(e.params)).z_c;
    const HopSizing size = plate.hopSizing(cascadeServoTau(e.params));

    auto accelAt = [&](double z) {
        return RollingBallDynamics::rolling_factor() * g *
               std::sin(tk.max_conditioned_tilt(z, kRatesUntrustworthyAbove));
    };

    // 15.63 deg and 1.888 m/s^2 at the home heave; 14.52 deg and 1.757 a hop's
    // charge below it.  **Seven per cent of the ball acceleration the plate can
    // command, which is what a hop gives up while it is charging.**
    ASSERT_NEAR(tk.max_conditioned_tilt(home_z, kRatesUntrustworthyAbove) * 180 / M_PI,
                15.63, 0.05);
    ASSERT_NEAR(accelAt(home_z), 1.888, 0.005);
    ASSERT_NEAR(tk.max_conditioned_tilt(home_z - size.drop_m,
                                        kRatesUntrustworthyAbove) * 180 / M_PI,
                14.52, 0.05);
    ASSERT_NEAR(accelAt(home_z - size.drop_m), 1.757, 0.005);
    ASSERT_TRUE(accelAt(home_z - size.drop_m) > 0.92 * accelAt(home_z));

    // And the plate really does stay inside that excursion — this is what
    // `HopCycle::z_arm` is for.  The loop moves the heave too, so the bound is
    // on the hop's share: twice the drop covers both.
    const Hop r = follow(e, defaultGain(e), openingCircle(), true, 30.0);
    ASSERT_TRUE(home_z - r.lowest_z < 2.0 * size.drop_m);
}

// **`bounceFloorSpeed` against a deliberate hop**, which is the interaction
// #23 left untested because nothing hopped on purpose yet.
//
// It does not terminate a commanded hop.  The margin is far smaller than the
// floor's own header implies, though, and the header's reasoning is what is
// wrong rather than its number:
//
//   - The floor bounds the RELATIVE rebound at an impact.  Its header compares
//     that to a hop's HEIGHT — 0.34 mm at 60 Hz against a 6.6 mm passive hop —
//     and the two are not the same quantity.  At a release the ball and the
//     plate are nearly comoving by construction, because the ball leaves
//     carrying the plate's velocity, so the relative speed the floor tests is
//     small at exactly the instant a hop begins even though the hop is seventy
//     times the height the floor can represent.
//   - And a throw DOES resolve an impact on the frame it lets go: the plate is
//     still accelerating up the ramp, so it runs into the ball it has just
//     released.  Measured, 7 of 8 throws at 60 Hz.  The floor is therefore
//     asked about every hop this controller commands.
//
// What it clears the floor by is a constant, and that is the structural part.
// The plate's rise gains `rise_rate dt / stroke_s = margin g dt / 2` a frame
// while the ramp runs, and the ball falls `g dt / 2` against it, so the
// approach is about `(1 + margin) g dt / 2` and the rebound `e` of that —
// against a floor of `g dt / 2`.  **The ratio is `e (1 + margin)`, about 2.4,
// and `dt` cancels**, which is why a finer frame does not widen it.  Measured:
// 1.8 at 60 Hz, 1.2 at 120 and 1.3 at 240, the shortfall being the heave solve
// yielding to the retreat and the ramp's first frame being partial.
//
// The edge was found by moving the floor rather than argued for: at twice
// `g dt / 2` the controller stops hopping at all.  So the margin is a factor of
// about two, it does not grow with the frame rate, and the only lever on it is
// `kHopRiseMargin` — which is bounded above at 1.7 by the carry-versus-strike
// step.  Recorded rather than widened.
void test_the_bounce_floor_cannot_end_a_commanded_hop() {
    const auto models = getBuiltinModels();
    const ModelEntry& e = cascadeModel(models);
    const SimPlate plate = cascadePlate(e.params);
    const double g = plate.gravity();
    const double tau = cascadeServoTau(e.params);
    const HopSizing size = plate.hopSizing(tau);

    double previous_ratio = 0.0;
    for (double hz : {60.0, 120.0, 240.0}) {
        const double dt = 1.0 / hz;
        const double floor = bounceFloorSpeed(g, dt);

        // The rate the throw commands is far above the floor at every frame
        // rate the hop works at, and that margin does grow as the frame
        // shrinks — it is the RELATIVE rebound below that does not.
        ASSERT_TRUE(size.rise_rate > 8.0 * floor);

        const double floor_hop = floor * floor / (2.0 * g);
        const double hop = size.rise_rate * size.rise_rate / (2.0 * g);
        const double ratio = hop / floor_hop;
        ASSERT_NEAR(ratio, 4.0 * kHopRiseMargin * kHopRiseMargin *
                               (tau / dt) * (tau / dt), 1e-6 * ratio);
        ASSERT_TRUE(ratio > previous_ratio);   // finer frames, wider margin
        previous_ratio = ratio;

        const Hop r = follow(e, defaultGain(e), openingCircle(), true, 30.0, dt);
        ASSERT_TRUE(r.thrown >= 6);

        // A throw DOES resolve an impact on the frame it lets go — the plate
        // is still climbing the ramp and runs into the ball it just released —
        // so the floor is asked about every hop this controller commands.  It
        // cleared it at every frame rate.
        ASSERT_TRUE(r.impacts_at_release > 0);
        ASSERT_TRUE(r.weakest_release < 1e8);   // so a release impact was seen
        ASSERT_TRUE(r.weakest_release > 1.15 * floor);
        // ...and by a factor that does NOT improve with the frame rate, which
        // is the point.  `e (1 + margin)` is 2.4 and the yielding takes it to
        // between 1.2 and 1.8; nothing here is allowed to drift above it
        // silently, because that would mean the release had turned back into a
        // strike.
        ASSERT_TRUE(r.weakest_release < kRestitution * (1.0 + kHopRiseMargin) * floor);

        // The hop itself clears the height the floor can end by two orders.
        ASSERT_TRUE(r.apex > 20.0 * floor_hop);
    }
}

// The flight is still the constraint's, and the constraint still holds.  The
// hop is what makes `p` positive at a contact; it must not make it positive
// during a flight, which is where `u_rebound = e u + p (1 + e)` would pump.
void test_a_commanded_hop_does_not_pump_its_own_train() {
    const auto models = getBuiltinModels();
    const ModelEntry& e = cascadeModel(models);
    for (const char* nm : {"Detuned", "Nominal", "Aggressive"}) {
        const Hop r = follow(e, gainForPreset(e, presetNamed(nm)),
                             openingCircle(), true, 30.0);
        ASSERT_TRUE(r.thrown >= 6);
        // Measured at 0.0002 m/s under every preset: the constraint holding
        // rather than yielding.  `holdContactDown` documents Aggressive giving
        // up as much as 0.296 under a `kMaxNudgeSpeed` shove, and a commanded
        // hop asks nothing like that of the servos.
        ASSERT_TRUE(r.worst_rise < 0.005);
        ASSERT_TRUE(!r.lost);
        ASSERT_TRUE(!r.changed_assembly);
        // Every preset hops, and to the same height: the hop is sized by the
        // plant rather than by the gain, which is what makes it the same
        // feature under all three.
        ASSERT_NEAR(r.apex, 0.0247, 0.006);
    }
}

// **The ball is never lost at any setting the sliders offer**, which is the
// same promise `test_trajectory` makes for tracking and the one thing a demo
// control may not break.  And **the ball still follows the path while hopping,
// to a tracking error measured at every one of them** — #23 asked for a stated
// error and a number that lives only in prose is the failure the decision
// record's D15 names.
//
// Nine of the eighteen never hop, and that is the feature working: a ball
// crossing the plate faster than `max_ball_speed` is one whose flight would not
// come down on the plate, so the throw is refused and the tracking is exactly
// what it was.  See `stepHop`.
//
// **The cornered shapes are in scope**, and this is where that is decided: the
// square and the triangle hop at their slow laps exactly as the circle does,
// and are refused at their fast ones exactly as the circle is.  The fillet
// (#31) is what makes that true — a reference with a step in its velocity at
// every corner would have the loop against its stops there, which is the second
// gate.
void test_no_offered_setting_loses_the_ball_to_a_hop() {
    const auto models = getBuiltinModels();
    const ModelEntry& e = cascadeModel(models);
    const SimPlate plate = cascadePlate(e.params);
    const Eigen::MatrixXd K = defaultGain(e);

    int hopped = 0, refused = 0;
    double worst_apex = 0.0, worst_rise = 0.0;
    double worst_hopping_error = 0.0, best_hopping_error = 1.0;
    double worst_cost = -1.0;   ///< the most a hop added to the tracking error
    for (PathShape shape : {PathShape::Circle, PathShape::Square, PathShape::Triangle}) {
        for (double radius : {0.06, 0.12, 0.18}) {
            for (int fast = 0; fast < 2; ++fast) {
                SetpointPath q;
                q.shape = shape;
                q.radius_m = radius;
                q.period_s = fast ? clampPeriod(plate.feasible(q), 0.0) : 30.0;

                const double seconds = 3.0 * q.period_s + 20.0;
                const Hop flat = follow(e, K, q, false, seconds);
                const Hop r = follow(e, K, q, true, seconds);
                ASSERT_TRUE(!r.lost);
                ASSERT_TRUE(!r.changed_assembly);

                if (r.thrown == 0) {
                    ++refused;
                    // A refused hop is not a degraded one: the command is
                    // untouched, so the tracking is the same run.
                    ASSERT_NEAR(r.mean_err, flat.mean_err, 1e-12);
                    continue;
                }
                ++hopped;
                worst_apex = std::max(worst_apex, r.apex);
                worst_rise = std::max(worst_rise, r.worst_rise);
                worst_hopping_error = std::max(worst_hopping_error, r.mean_err);
                best_hopping_error = std::min(best_hopping_error, r.mean_err);
                worst_cost = std::max(worst_cost, r.mean_err - flat.mean_err);
            }
        }
    }
    ASSERT_EQ(hopped, 9);
    ASSERT_EQ(refused, 9);

    // The hop is the same size wherever it happens, and never pumps.
    ASSERT_TRUE(worst_apex < 0.030);
    ASSERT_TRUE(worst_rise < 0.005);

    // **The stated tracking error, over every setting that hops**: 1.8 to
    // 3.9 mm, against 3.7 to 4.5 mm for the same nine runs without the hop.
    // Hopping costs nothing measurable here and on most of them it is better —
    // the tilt is untouched by construction, and `predictedLanding` aiming the
    // plate at where the ball will come down cancels a lag the rolling loop
    // integrates.  The bound below is on the COST rather than on the error, so
    // a path that is simply hard to track cannot be mistaken for a hop that
    // made it harder.
    ASSERT_TRUE(best_hopping_error > 0.0015);
    ASSERT_TRUE(worst_hopping_error < 0.0040);
    ASSERT_TRUE(worst_cost < 0.0005);
}

// **The disturbances are live alongside the hop**, and the two together are
// what a visitor will actually do: tick the box, then shove the ball.  The
// sweep above is undisturbed laps, so on its own it says nothing about that.
//
// `kMaxNudgeSpeed` is the hardest shove the interface offers, and #23 measured
// it separating the ball on about one run in nine of a lap-by-direction grid
// with no hop at all.  With the hop on, the shove and the throw can land
// together — and they must still not lose the ball.
void test_a_shove_during_a_hop_still_keeps_the_ball() {
    const auto models = getBuiltinModels();
    const ModelEntry& e = cascadeModel(models);
    const SimPlate plate = cascadePlate(e.params);
    const Eigen::MatrixXd K = defaultGain(e);

    SimInput in;
    in.design = cascadeDesign(e.params);
    in.design.K = K;
    in.path = plate.feasible(openingCircle());
    in.hop_enabled = true;

    int shoves = 0, lost = 0, assembly = 0, hopped = 0;
    double worst_radius = 0.0;
    // Eight points of the lap crossed with eight directions, the shove timed
    // to fall at a different point of the hop cycle each round.
    for (int lap_point = 0; lap_point < 8; ++lap_point) {
        for (int dir = 0; dir < 8; ++dir) {
            const Eigen::Vector2d start = pathPoint(in.path, 0.0);
            SimState s = simStart(plate, in.design.home_leg_rad,
                                  Eigen::Vector4d(start(0), start(1), 0.0, 0.0));
            const double shove_at =
                in.path.period_s * (1.0 + lap_point / 8.0) + dir * in.dt;
            const double theta = dir * 2.0 * M_PI / 8.0;

            bool shoved = false;
            double t = 0.0;
            for (int k = 0; k < static_cast<int>(35.0 / in.dt); ++k) {
                if (!shoved && t >= shove_at && !s.ball.airborne) {
                    s.ball.rolling(2) += kMaxNudgeSpeed * std::cos(theta);
                    s.ball.rolling(3) += kMaxNudgeSpeed * std::sin(theta);
                    shoved = true;
                    ++shoves;
                }
                const SimReport f = stepSim(plate, in, s);
                t += in.dt;
                worst_radius = std::max(
                    worst_radius, std::hypot(f.ball_plate(0), f.ball_plate(1)));
                if (f.left_plate) { ++lost; break; }
                if (!onBuiltAssembly(plate.kinematics(), s.alpha_rad, s.pose))
                    ++assembly;
            }
            if (s.hop.thrown > 0) ++hopped;
        }
    }

    ASSERT_EQ(shoves, 64);
    // Measured: none of the 64 loses the ball, none changes assembly, every
    // one of them hops, and the ball reaches 164 mm against a 280 mm rim.  The
    // shove is what takes it out there — 120 mm of path plus a 0.20 m/s kick —
    // and a shove arriving while the ball is being thrown or is mid-train is
    // the interaction the undisturbed sweep cannot reach at all.
    ASSERT_EQ(lost, 0);
    ASSERT_EQ(assembly, 0);
    ASSERT_TRUE(worst_radius < 0.18);
    ASSERT_EQ(hopped, 64);
}

}  // namespace

int main() {
    test_the_hop_is_sized_by_the_servo_the_ball_and_the_rim();
    test_disabling_the_hop_forgets_the_cycle();
    test_a_hop_charges_down_before_it_throws_up();
    test_the_strokes_ramp_rather_than_step();
    test_the_strokes_end_on_travel_with_the_timer_as_a_backstop();
    test_a_failed_throw_does_not_march_the_plate_down();
    test_the_flight_belongs_to_the_constraint();
    test_only_a_throw_that_let_go_is_counted();
    test_a_throw_is_refused_while_the_ball_is_fast_or_the_loop_is_on_its_stops();

    test_the_heave_solve_hits_its_target_at_either_sign();
    test_a_throw_moves_all_three_legs_the_same_way();
    test_the_heave_solve_declines_what_the_constraint_declines();

    test_the_plate_hops_the_ball_while_it_tracks();
    test_the_hop_spends_a_stated_amount_of_tilt_authority();
    test_the_bounce_floor_cannot_end_a_commanded_hop();
    test_a_commanded_hop_does_not_pump_its_own_train();
    test_no_offered_setting_loses_the_ball_to_a_hop();
    test_a_shove_during_a_hop_still_keeps_the_ball();

    std::printf("test_hop_drive: all passed\n");
    return 0;
}
