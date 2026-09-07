// tests/test_auto_balance.cpp
//
// Closing the loop is the one place where the two halves of this application
// have to agree numerically rather than merely coexist.  The gain comes from a
// LINEAR model; it drives the NONLINEAR simulator — a Newton-solved mechanism
// and a rolling ball, neither of which the design ever saw.
//
// So the tests come in two layers.  The first pins the wiring: which state
// goes in which slot, what the setpoint means, where the servo travel clips.
// The second runs the real thing and asserts what the ticket actually asks
// for — that the ball balances under a sensible tuning and visibly does not
// under a poor one.  A criterion phrased as "visibly" is still a measurable
// claim, and measuring it here is what keeps it from being re-judged by eye
// every time the model moves.
#include "auto_balance.h"

#include "ball_contact.h"
#include "ball_sim.h"
#include "cascade_fixture.h"
#include "table_kinematics.h"
#include "test_helpers.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>

using namespace caliburn;

namespace {

constexpr double kDeg = M_PI / 180.0;
constexpr double kHome = 45.0 * kDeg;
constexpr double kBallRadius = kFixtureBallRadius;
constexpr double kBallFriction = kFixtureBallFriction;

AutoBalanceDesign designWith(const Eigen::MatrixXd& K) {
    AutoBalanceDesign d;
    d.K = K;
    d.home_leg_rad = kHome;
    d.servo_tau = 0.05;
    return d;
}

// ---------------------------------------------------------------------------
// Layer 1: the wiring
// ---------------------------------------------------------------------------

void test_state_vector_layout() {
    const std::array<double, 3> alpha = {kHome + 0.1, kHome - 0.2, kHome};
    const Eigen::Vector4d ball(0.03, -0.04, 0.5, -0.6);

    const Eigen::VectorXd x = cascadeState(alpha, kHome, ball);

    ASSERT_EQ((int)x.size(), 7);
    ASSERT_NEAR(x(0), 0.1, 1e-12);
    ASSERT_NEAR(x(1), -0.2, 1e-12);
    ASSERT_NEAR(x(2), 0.0, 1e-12);
    ASSERT_NEAR(x(3), 0.03, 1e-12);   // ball x
    ASSERT_NEAR(x(4), -0.04, 1e-12);  // ball y
    ASSERT_NEAR(x(5), 0.5, 1e-12);    // ball x'
    ASSERT_NEAR(x(6), -0.6, 1e-12);   // ball y'
}

// The deviation is from the DESIGN's home angle, not from any fixed 45 degrees.
// Get this wrong and the loop still converges — to a tilted plate.
void test_state_deviation_follows_the_design_home() {
    const std::array<double, 3> alpha = {60.0 * kDeg, 60.0 * kDeg, 60.0 * kDeg};
    const Eigen::VectorXd x = cascadeState(alpha, 60.0 * kDeg, Eigen::Vector4d::Zero());
    for (int i = 0; i < 3; ++i) ASSERT_NEAR(x(i), 0.0, 1e-12);
}

void test_shape_gate() {
    AutoBalanceDesign d = designWith(Eigen::MatrixXd::Zero(3, 7));
    ASSERT_TRUE(gainFitsCascade(d));
    d.K = Eigen::MatrixXd::Zero(1, 4);
    ASSERT_TRUE(!gainFitsCascade(d));
    d.K = Eigen::MatrixXd();
    ASSERT_TRUE(!gainFitsCascade(d));
}

// A ball at rest at the setpoint, legs at home: nothing to do.  This is the
// claim that makes a position setpoint need no feedforward.
// The model panel carries geometry as `float` sliders; the plate holds
// `double` literals.  0.300f widened is 0.3000000119, so a comparison tight
// enough to be called exact refuses two mechanisms that are the same object —
// and the balance loop then never becomes available at all.
void test_plant_identity_survives_the_float_round_trip() {
    TableParams plate;
    plate.R_ground = 0.300;
    plate.R_table = 0.300;
    plate.L1 = 0.150;
    plate.L2 = 0.150;

    // Exactly how the model panel's float sliders reach `cascadeMechanism`.
    const auto models = getBuiltinModels();
    const TableParams from_sliders = cascadeMechanism(cascadeModel(models).params);
    ASSERT_TRUE(samePlant(plate, 9.81, from_sliders, 9.81));

    // A millimetre of leg, however, is a different plate.
    TableParams longer = from_sliders;
    longer.L1 = 0.151;
    ASSERT_TRUE(!samePlant(plate, 9.81, longer, 9.81));

    // And so is the same plate on the moon — the gain would be solved for an
    // acceleration the simulated ball never feels.
    ASSERT_TRUE(!samePlant(plate, 9.81, from_sliders, 1.62));

    // A design nobody filled in is all zeros, and matches nothing.
    ASSERT_TRUE(!samePlant(plate, 9.81, TableParams{}, 0.0));
}

void test_zero_error_commands_the_home_pose() {
    Eigen::MatrixXd K = Eigen::MatrixXd::Random(3, 7) * 10.0;
    const AutoBalanceDesign d = designWith(K);

    const std::array<double, 3> alpha = {kHome, kHome, kHome};
    const Eigen::Vector4d ball(0.07, -0.02, 0.0, 0.0);

    const LegCommand c = legCommand(cascadeKinematics(), d, alpha, ball, {{0.07, -0.02}});
    for (int i = 0; i < 3; ++i) ASSERT_NEAR(c.alpha_rad[i], kHome, 1e-12);
    ASSERT_TRUE(!c.saturated);
}

void test_command_is_home_minus_k_error() {
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(3, 7);
    K(0, 3) = 2.0;   // leg 0 reacts to ball x
    K(1, 4) = -3.0;  // leg 1 reacts to ball y
    K(2, 0) = 0.5;   // leg 2 reacts to its own neighbour's deviation
    const AutoBalanceDesign d = designWith(K);

    const std::array<double, 3> alpha = {kHome + 0.04, kHome, kHome};
    const Eigen::Vector4d ball(0.01, 0.02, 0.0, 0.0);

    const LegCommand c = legCommand(cascadeKinematics(), d, alpha, ball, {});
    ASSERT_NEAR(c.alpha_rad[0], kHome - 2.0 * 0.01, 1e-12);
    ASSERT_NEAR(c.alpha_rad[1], kHome + 3.0 * 0.02, 1e-12);
    ASSERT_NEAR(c.alpha_rad[2], kHome - 0.5 * 0.04, 1e-12);
    ASSERT_TRUE(!c.saturated);
}

// The setpoint enters as a reference STATE, so displacing both the ball and the
// setpoint together is a no-op while displacing only one is not.
void test_setpoint_moves_the_target() {
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(3, 7);
    K(0, 3) = 1.0;
    const AutoBalanceDesign d = designWith(K);
    const std::array<double, 3> alpha = {kHome, kHome, kHome};

    const LegCommand at = legCommand(cascadeKinematics(), d, alpha,
                                     Eigen::Vector4d(0.05, 0, 0, 0), {{0.05, 0.0}});
    ASSERT_NEAR(at.alpha_rad[0], kHome, 1e-12);

    const LegCommand off = legCommand(cascadeKinematics(), d, alpha,
                                      Eigen::Vector4d(0.05, 0, 0, 0), {});
    ASSERT_NEAR(off.alpha_rad[0], kHome - 0.05, 1e-12);
}

// The reference VELOCITY, which is the other half of the reference state and
// the whole of the trajectory feedforward.  It enters exactly as the position
// does — as a term of `x_ref`, differenced against the measurement — so a ball
// moving at the speed the setpoint is moving at is a ball with no error to
// answer, and the plate stays level.
//
// It lives here, in the loop, rather than in whoever is driving the setpoint.
// It used to be applied by biasing the measurement handed in, two files away
// in `plate_view`, with a third copy in `test_trajectory`; the arithmetic was
// identical and the arrangement was not, because the loop's own header went on
// saying it needed no feedforward.  See #24.
void test_the_reference_velocity_enters_as_a_velocity_error() {
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(3, 7);
    K(0, 5) = 4.0;   // leg 0 reacts to the ball's x velocity
    K(1, 6) = -2.0;  // leg 1 to its y velocity
    const AutoBalanceDesign d = designWith(K);
    const std::array<double, 3> alpha = {kHome, kHome, kHome};

    // Ball running at the setpoint's own speed: no error, level plate.
    BallReference ref;
    ref.velocity = Eigen::Vector2d(0.08, -0.03);
    const LegCommand with = legCommand(cascadeKinematics(), d, alpha,
                                       Eigen::Vector4d(0, 0, 0.08, -0.03), ref);
    for (int i = 0; i < 3; ++i) ASSERT_NEAR(with.alpha_rad[i], kHome, 1e-12);

    // The same ball with a HELD setpoint is a ball the loop tries to stop.
    const LegCommand without = legCommand(cascadeKinematics(), d, alpha,
                                          Eigen::Vector4d(0, 0, 0.08, -0.03), {});
    ASSERT_NEAR(without.alpha_rad[0], kHome - 4.0 * 0.08, 1e-12);
    ASSERT_NEAR(without.alpha_rad[1], kHome + 2.0 * -0.03, 1e-12);
}

// What the tuning panel says about the feedforward, as arithmetic.
//
// The feedforward has no gain of its own: the reference velocity is multiplied
// by K's velocity columns, and LQR designs those from the `x' ball` and
// `y' ball` weights.  So those two sliders DO set how hard the loop chases a
// moving setpoint, and the panel says so — which is a claim, and this is the
// test under it.  Inventing a separate feedforward gain would be a second
// opinion about a number K already owns.
void test_the_ball_velocity_weights_set_the_feedforwards_strength() {
    const auto models = getBuiltinModels();
    const auto& e = cascadeModel(models);

    Eigen::VectorXd q = defaultLqrStateWeights(7);
    const Eigen::VectorXd r = defaultLqrInputWeights(3);
    const Eigen::MatrixXd soft = gainFor(e, q, r);

    q(5) = q(6) = 10.0 * q(5);          // the two sliders the panel names
    const Eigen::MatrixXd hard = gainFor(e, q, r);

    // Only the velocity columns are what the feedforward is multiplied by.
    ASSERT_TRUE(hard.col(5).norm() > 2.0 * soft.col(5).norm());
    ASSERT_TRUE(hard.col(6).norm() > 2.0 * soft.col(6).norm());

    // And so the same moving setpoint asks the legs for more.
    BallReference ref;
    ref.velocity = Eigen::Vector2d(0.075, 0.0);   // the opening circle's speed
    const std::array<double, 3> alpha = {kHome, kHome, kHome};
    const Eigen::Vector4d still(0, 0, 0, 0);      // a ball not yet keeping up

    AutoBalanceDesign ds = designWith(soft), dh = designWith(hard);
    const LegCommand cs = legCommand(cascadeKinematics(), ds, alpha, still, ref);
    const LegCommand ch = legCommand(cascadeKinematics(), dh, alpha, still, ref);

    double span_s = 0.0, span_h = 0.0;
    for (int i = 0; i < 3; ++i) {
        span_s = std::max(span_s, std::abs(cs.alpha_rad[i] - kHome));
        span_h = std::max(span_h, std::abs(ch.alpha_rad[i] - kHome));
    }
    ASSERT_TRUE(span_h > span_s);
}

// Travel limits clamp; the workspace then clips.  Two stages, because they
// are two different constraints: a leg has a stop, and a triple of legs has to
// describe a plate that exists.  Since #22 the second one is enforced too, so
// a command driven at both stops does not come back sitting on them.
void test_command_clamps_to_servo_travel_then_to_the_workspace() {
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(3, 7);
    K(0, 3) = 1e4;
    K(1, 3) = -1e4;
    const AutoBalanceDesign d = designWith(K);
    const std::array<double, 3> alpha = {kHome, kHome, kHome};

    const LegCommand c = legCommand(cascadeKinematics(), d, alpha,
                                    Eigen::Vector4d(0.05, 0, 0, 0), {});
    ASSERT_TRUE(c.saturated);
    ASSERT_TRUE(c.clipped_to_workspace);

    // Both stops asked for a 70 degree spread, which no assembly of this
    // mechanism has.  What comes back is inside the travel rather than on it.
    ASSERT_TRUE(c.alpha_rad[0] > d.alpha_min_rad);
    ASSERT_TRUE(c.alpha_rad[1] < d.alpha_max_rad);

    // Direction survives: leg 0 still driven down, leg 1 up, leg 2 untouched.
    ASSERT_TRUE(c.alpha_rad[0] < kHome);
    ASSERT_TRUE(c.alpha_rad[1] > kHome);
    ASSERT_NEAR(c.alpha_rad[2], kHome, 1e-12);

    // And the plate can be built at what it asked for.
    const TableKinematics& tk = cascadeKinematics();
    const double mean = (c.alpha_rad[0] + c.alpha_rad[1] + c.alpha_rad[2]) / 3.0;
    ASSERT_TRUE(tk.solve_pose(c.alpha_rad, tk.home_pose(mean)).converged);
}

// The clamp on its own, where the clamped triple IS assemblable: the limit is
// still a hard stop, and nothing gets scaled back for no reason.
//
// The travel has to be narrowed for this case to exist at all, which is worth
// knowing on its own.  The shipped limits are 10 and 80 degrees, and a SINGLE
// leg taken to either of them with the other two at home already has no
// assembly — the mechanism's real range about the home pose is far narrower
// than its servos'.  That is why the workspace clip is doing work on almost
// every saturated frame rather than in a rare corner.
void test_a_reachable_clamp_is_left_on_the_stop() {
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(3, 7);
    K(0, 3) = 1e4;                       // one leg only, driven to its stop
    AutoBalanceDesign d = designWith(K);
    d.alpha_min_rad = 40.0 * kDeg;
    d.alpha_max_rad = 50.0 * kDeg;
    const std::array<double, 3> alpha = {kHome, kHome, kHome};

    const LegCommand c = legCommand(cascadeKinematics(), d, alpha,
                                    Eigen::Vector4d(0.05, 0, 0, 0), {});
    ASSERT_TRUE(c.saturated);
    ASSERT_TRUE(!c.clipped_to_workspace);
    ASSERT_NEAR(c.alpha_rad[0], d.alpha_min_rad, 1e-12);
    ASSERT_NEAR(c.alpha_rad[1], kHome, 1e-12);
    ASSERT_NEAR(c.alpha_rad[2], kHome, 1e-12);
}

void test_wrong_shape_commands_the_home_pose() {
    const AutoBalanceDesign d = designWith(Eigen::MatrixXd::Ones(2, 3));
    const std::array<double, 3> alpha = {kHome + 0.2, kHome, kHome};
    const LegCommand c = legCommand(cascadeKinematics(), d, alpha,
                                    Eigen::Vector4d(0.1, 0.1, 1, 1), {});
    for (int i = 0; i < 3; ++i) ASSERT_NEAR(c.alpha_rad[i], kHome, 1e-12);
    ASSERT_TRUE(!c.saturated);
}

void test_servo_lag_is_first_order() {
    const double tau = 0.05;
    const std::array<double, 3> cmd = {1.0, 1.0, 1.0};
    std::array<double, 3> a = {0.0, 0.0, 0.0};

    // One time constant, reached in one step and in many: the exact form makes
    // the answer independent of how the interval was subdivided.
    const std::array<double, 3> one_shot = stepServos(a, cmd, tau, tau);
    ASSERT_NEAR(one_shot[0], 1.0 - std::exp(-1.0), 1e-12);

    for (int k = 0; k < 100; ++k) a = stepServos(a, cmd, tau, tau / 100.0);
    ASSERT_NEAR(a[0], 1.0 - std::exp(-1.0), 1e-9);

    // No step, no motion — and no division by a zero tau.
    const std::array<double, 3> still = stepServos(a, cmd, tau, 0.0);
    ASSERT_NEAR(still[0], a[0], 1e-15);
    const std::array<double, 3> degenerate = stepServos(a, cmd, 0.0, 1.0 / 60.0);
    ASSERT_TRUE(std::isfinite(degenerate[0]));

    // Each leg keeps its own command.
    const std::array<double, 3> mixed =
        stepServos({0.0, 0.0, 0.0}, {1.0, -1.0, 0.0}, tau, tau);
    ASSERT_TRUE(mixed[0] > 0.0 && mixed[1] < 0.0);
    ASSERT_NEAR(mixed[2], 0.0, 1e-15);
}

// ---------------------------------------------------------------------------
// Layer 2: the loop, closed around the nonlinear plate
// ---------------------------------------------------------------------------

struct SimResult {
    double final_radius;   // [m] distance of the ball from its setpoint
    double peak_radius;    // [m]
    double settle_time;    // [s] after which it never again left 5 mm; -1 = never
    double final_tilt_deg; // [deg] the larger of |roll| and |pitch| at the end
    bool saturated;        // a leg command hit a travel limit at least once
    bool clipped;          // ...or asked for a pose the mechanism cannot make
    bool left_plate;
};

// The same four calls, in the same order, that PlateView::step makes:
// legCommand -> stepServosOnPlate -> solve_pose -> stepBall.  Nothing here is
// linearised, and nothing here knows what Q and R were.
SimResult runClosedLoop(const ModelEntry& e,
                        const Eigen::MatrixXd& K,
                        const Eigen::Vector4d& ball0,
                        double x_sp,
                        double y_sp,
                        double duration) {
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
    ball.rolling = ball0;
    TablePose pose = tk.home_pose(d.home_leg_rad);
    PlateMotion pm = plateMotion(tk, pose, alpha, {0.0, 0.0, 0.0}), pm_prev = pm;
    const double g = cascadeGravity(e.params);

    const double dt = 1.0 / 60.0;
    const int steps = static_cast<int>(duration / dt);

    SimResult r{0.0, 0.0, -1.0, 0.0, false, false, false};
    for (int k = 0; k < steps; ++k) {
        const Eigen::Matrix<double, 6, 1> bp = plateFrame(ball, pm, kBallRadius);
        Eigen::Vector4d seen(bp(0), bp(1), bp(3), bp(4));
        if (ball.airborne) {
            const Eigen::Vector2d land = predictedLanding(bp, kBallRadius, g);
            seen << land(0), land(1), 0.0, 0.0;
        }
        const LegCommand c = legCommand(tk, d, alpha, seen, {{x_sp, y_sp}});
        if (c.saturated) r.saturated = true;
        if (c.clipped_to_workspace) r.clipped = true;
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
        const Eigen::Matrix<double, 6, 1> now = plateFrame(ball, pm, kBallRadius);

        const double radius = std::hypot(now(0) - x_sp, now(1) - y_sp);
        r.peak_radius = std::max(r.peak_radius, radius);
        if (radius > 0.005) r.settle_time = (k + 1) * dt;
        if (!ballOnPlate(Eigen::Vector4d(now(0), now(1), now(3), now(4)),
                         tp.R_table, kBallRadius)) {
            r.left_plate = true;
            r.final_radius = radius;
            return r;
        }
    }
    const Eigen::Matrix<double, 6, 1> fin = plateFrame(ball, pm, kBallRadius);
    r.final_radius = std::hypot(fin(0) - x_sp, fin(1) - y_sp);
    r.final_tilt_deg = std::max(std::abs(pose.phi), std::abs(pose.theta)) / kDeg;
    if (r.settle_time >= duration - dt) r.settle_time = -1.0;
    return r;
}

// 60 mm out along x and 40 mm back along y, at rest: a quarter of the plate's
// radius, plainly visible in the 3D view and well outside anything a
// linearisation could excuse.
//
// The application used to OPEN on this displacement, which is why it was the
// one measured here.  It no longer does — that step input was the opening
// slam, and the demo now starts the ball already tracing a circle (see
// `attract_mode.h`).  The displacement stays as the yardstick anyway: "recover
// from a quarter of the plate, at rest" is the question this file exists to
// ask of a gain, and it is a harder one than the opening now poses.
const Eigen::Vector4d kDisplaced(0.06, -0.04, 0.0, 0.0);

// The weights the application opens on.  If THIS does not balance the ball,
// nobody who switches the controller to LQR ever sees the product work.
void test_default_weights_balance_the_ball() {
    const auto models = getBuiltinModels();
    const auto& e = cascadeModel(models);

    const SimResult s = runClosedLoop(e, defaultGain(e), kDisplaced, 0.0, 0.0, 20.0);

    ASSERT_TRUE(!s.left_plate);
    ASSERT_TRUE(s.settle_time > 0.0);
    ASSERT_TRUE(s.settle_time < 5.0);
    ASSERT_TRUE(s.final_radius < 0.005);
}

// The nudge buttons in Plate Control, as a number: an impulse in velocity,
// rejected.  This is the disturbance the "visibly degrades" criterion is read
// against, so the good case has to be pinned first.
void test_default_weights_reject_a_nudge() {
    const auto models = getBuiltinModels();
    const auto& e = cascadeModel(models);

    const SimResult s = runClosedLoop(
        e, defaultGain(e), Eigen::Vector4d(0.0, 0.0, 0.35, 0.25), 0.0, 0.0, 20.0);

    ASSERT_TRUE(!s.left_plate);
    ASSERT_TRUE(s.peak_radius > 0.05);      // it really was kicked
    ASSERT_TRUE(s.settle_time > 0.0);
    ASSERT_TRUE(s.settle_time < 5.0);
    ASSERT_TRUE(s.final_radius < 0.005);
}

// A ball at rest anywhere on a flat plate is an equilibrium, so a position
// setpoint is reachable with no feedforward at all.  That claim is arithmetic
// in `legCommand`; here it is the nonlinear plate that has to agree.
void test_setpoint_is_tracked() {
    const auto models = getBuiltinModels();
    const auto& e = cascadeModel(models);

    const SimResult s = runClosedLoop(
        e, defaultGain(e), Eigen::Vector4d::Zero(), 0.08, 0.05, 20.0);

    ASSERT_TRUE(!s.left_plate);
    ASSERT_TRUE(s.settle_time > 0.0);
    ASSERT_TRUE(s.settle_time < 5.0);
    ASSERT_TRUE(s.final_radius < 0.005);
}

// "Moving the poles toward the imaginary axis should visibly make the ball
// sluggish" — the ticket's own words, as a comparison the user can reproduce
// by dragging the ball-position weights down.
void test_looser_weights_are_visibly_sluggish() {
    const auto models = getBuiltinModels();
    const auto& e = cascadeModel(models);

    Eigen::VectorXd q_loose(7);
    q_loose << 1, 1, 1, 50, 50, 5, 5;
    const Eigen::MatrixXd K_loose = gainFor(e, q_loose, defaultLqrInputWeights(3));

    const SimResult loose = runClosedLoop(e, K_loose, kDisplaced, 0.0, 0.0, 20.0);
    const SimResult base = runClosedLoop(e, defaultGain(e), kDisplaced, 0.0, 0.0, 20.0);

    ASSERT_TRUE(loose.settle_time > 0.0);
    ASSERT_TRUE(base.settle_time > 0.0);
    // Not "a bit slower": four times slower, which is the difference between
    // watching it arrive and wondering whether it is moving.
    ASSERT_TRUE(loose.settle_time > 4.0 * base.settle_time);
}

// The third acceptance criterion, as a number.  Legs priced a thousand times
// above the ball buys a controller that will barely move them: it stalls the
// plate inside the friction dead band and the ball simply stays where it was
// put, or wanders half the plate away when it is kicked.
//
// R = 1000 is the top of the model panel's own slider, deliberately.  A test
// that proved degradation at 1e5 would be proving it about a tuning nobody can
// reach through the UI, which is not the criterion.  Every number below is
// against the good case measured in the two tests above: 1.1 mm and 1.8 s from
// the same displacement, 81 mm of peak and 2.2 s from the same kick.
void test_poor_tuning_visibly_degrades() {
    const auto models = getBuiltinModels();
    const auto& e = cascadeModel(models);

    const Eigen::MatrixXd K_poor = gainFor(
        e, defaultLqrStateWeights(7), Eigen::VectorXd::Constant(3, 1000.0));

    const SimResult held = runClosedLoop(e, K_poor, kDisplaced, 0.0, 0.0, 20.0);
    ASSERT_TRUE(held.settle_time < 0.0);    // never settles, in twenty seconds
    ASSERT_TRUE(held.final_radius > 0.05);  // and never really tried
    // Stalled inside the dead band rather than converging slowly: the tilt it
    // is asking for is one the ball cannot be started by.
    ASSERT_TRUE(held.final_tilt_deg < std::atan(kBallFriction) / kDeg);

    const SimResult kicked = runClosedLoop(
        e, K_poor, Eigen::Vector4d(0.0, 0.0, 0.35, 0.25), 0.0, 0.0, 20.0);
    ASSERT_TRUE(kicked.peak_radius > 0.12);   // half the plate, against 81 mm
    ASSERT_TRUE(kicked.settle_time < 0.0);
    ASSERT_TRUE(kicked.final_radius > 0.05);
}

// Why the defaults above are not unit weights, kept as a fact rather than a
// memory.  The rolling model's Coulomb resistance is a dead band: below
// atan(c_rr) of tilt the ball cannot be started at all, and a state feedback
// with no integral term parks inside it.  A future reader will meet this as
// "the ball stops 37 mm off centre and the plate just sits there tilted", and
// it is a property of the plant, not a defect in the loop.
void test_unit_weights_park_in_the_friction_dead_band() {
    const auto models = getBuiltinModels();
    const auto& e = cascadeModel(models);

    const Eigen::MatrixXd K_unit =
        gainFor(e, Eigen::VectorXd::Ones(7), Eigen::VectorXd::Ones(3));
    const SimResult s = runClosedLoop(e, K_unit, kDisplaced, 0.0, 0.0, 20.0);

    ASSERT_TRUE(!s.left_plate);
    ASSERT_TRUE(s.settle_time < 0.0);        // stalls, well short of centre
    ASSERT_TRUE(s.final_radius > 0.02);
    // Stalled against the tilt that rolling resistance exactly cancels — to
    // within a tenth of it, because sign(v) leaves the ball creeping rather
    // than stopping dead.  An order of magnitude below the 8 deg the working
    // tuning uses, which is the whole difference between the two runs.
    ASSERT_NEAR(s.final_tilt_deg, std::atan(kBallFriction) / kDeg, 0.06);
}

// ---------------------------------------------------------------------------
// Layer 3: the preset tunings (#19)
// ---------------------------------------------------------------------------
//
// Three named tunings exist so that somebody who has never met a state
// weighting matrix can still watch what controller design does.  That makes
// every one of them a claim about the ball's behaviour, and the claims are
// printed on screen next to the buttons — so they are asserted here rather
// than re-judged by eye whenever the plant moves.

SimResult runPreset(const LqrPreset& p) {
    const auto models = getBuiltinModels();
    const auto& e = cascadeModel(models);
    return runClosedLoop(e, gainForPreset(e, p), kDisplaced, 0.0, 0.0, 30.0);
}

/// The Nominal run, computed once: it is the yardstick every other preset is
/// read against, and three tests want it.  `nominalPreset()` rather than
/// `defaultGain`, which is the same matrix by the invariant the test below
/// pins — going through the preset is what makes that a shared claim rather
/// than a coincidence this file relies on.
const SimResult& nominalRun() {
    static const SimResult r = runPreset(nominalPreset());
    return r;
}

// "Nominal is the tuning loaded at startup" — arithmetic, not a coincidence
// two sets of literals have to keep agreeing about.  `defaultLqrStateWeights`
// is DEFINED as the nominal preset, and this is the test that would notice if
// somebody inverted that dependency again.
void test_the_nominal_preset_is_the_startup_tuning() {
    const Eigen::VectorXd q = presetStateWeights(nominalPreset(), 7);
    const Eigen::VectorXd r = presetInputWeights(nominalPreset(), 3);
    ASSERT_TRUE(q.isApprox(defaultLqrStateWeights(7)));
    ASSERT_TRUE(r.isApprox(defaultLqrInputWeights(3)));

    // And it is the middle button, so a visitor arriving at the panel is
    // sitting on the preset the demo is already running — sluggish, shipped,
    // saturating, in that order.  Compared by address: `nominalPreset()` hands
    // back a reference INTO the table, and asserting that identity is what
    // stops a fourth preset being inserted above it unnoticed.
    ASSERT_TRUE(&lqrPresets()[1] == &nominalPreset());
}

// Only the ball's own weights vary between presets: the leg weights are not
// what the contrast is about, and a visitor who watched two knobs move could
// not say which one they had just seen the consequence of.
void test_the_presets_differ_only_in_the_ball_weights() {
    for (const auto& p : lqrPresets()) {
        const Eigen::VectorXd q = presetStateWeights(p, 7);
        ASSERT_NEAR(q(0), 1.0, 1e-12);
        ASSERT_NEAR(q(1), 1.0, 1e-12);
        ASSERT_NEAR(q(2), 1.0, 1e-12);
        ASSERT_NEAR(q(3), q(4), 1e-12);
        ASSERT_NEAR(q(5), q(6), 1e-12);
        ASSERT_NEAR(presetInputWeights(p, 3)(0), 1.0, 1e-12);
    }

    // A plant that is not the cascade gets unit weights, exactly as the
    // startup default does — the UI does not offer the presets for one, and
    // this is what stops a preset meaning something arbitrary if it ever did.
    ASSERT_TRUE(presetStateWeights(presetNamed("Aggressive"), 4)
                    .isApprox(Eigen::VectorXd::Ones(4)));
}

// Which button is lit.  The presets are an offer rather than a mode, so the
// answer is asked of the weights every frame rather than stored — and the
// interesting half is that dragging any slider puts it back to "none", which
// is the UI telling the truth about a tuning that is no longer a preset.
void test_the_lit_button_follows_the_weights() {
    const auto& presets = lqrPresets();
    for (int i = 0; i < static_cast<int>(presets.size()); ++i) {
        ASSERT_EQ(activePreset(presetStateWeights(presets[i], 7),
                               presetInputWeights(presets[i], 3)), i);
    }

    // Startup sits on Nominal, so a visitor arriving at the panel can see
    // which of the three they are already watching.
    ASSERT_EQ(activePreset(defaultLqrStateWeights(7), defaultLqrInputWeights(3)), 1);

    // One slider moved is no longer a preset, and neither is a plant that is
    // not the cascade — off it the three presets are the same unit weights,
    // and reporting a match would only be reporting that they stopped
    // differing.
    Eigen::VectorXd q = defaultLqrStateWeights(7);
    q(3) *= 1.1;
    ASSERT_EQ(activePreset(q, defaultLqrInputWeights(3)), -1);
    ASSERT_EQ(activePreset(defaultLqrStateWeights(7),
                           Eigen::VectorXd::Constant(3, 5.0)), -1);
    ASSERT_EQ(activePreset(Eigen::VectorXd::Ones(4), Eigen::VectorXd::Ones(2)), -1);
}

// The nominal run, as the yardstick the other two are read against.  Also the
// half of the aggressive claim that cannot be made by the aggressive run
// alone: "living against the servo stops" is only a contrast if the shipped
// tuning does not.
void test_the_nominal_preset_settles_without_touching_the_stops() {
    const SimResult s = nominalRun();
    ASSERT_TRUE(!s.left_plate);
    ASSERT_TRUE(!s.saturated);
    ASSERT_TRUE(s.settle_time > 1.5 && s.settle_time < 2.2);   // 1.80 measured
    ASSERT_TRUE(s.final_radius < 0.005);
}

// "Detuned is recognisably poor without being unstable" — the ticket's own
// words, and the two halves need separate assertions.  Poor: six times the
// settling time, which is the difference between watching it arrive and
// wondering whether it is moving.  Not unstable: it does arrive, to half a
// millimetre, and it never goes near a travel limit doing it.
void test_the_detuned_preset_is_sluggish_but_never_unstable() {
    const SimResult s = runPreset(presetNamed("Detuned"));
    const SimResult base = nominalRun();

    ASSERT_TRUE(!s.left_plate);
    ASSERT_TRUE(s.settle_time > 4.0 * base.settle_time);
    ASSERT_TRUE(s.final_radius < 0.005);
    ASSERT_TRUE(!s.saturated);
    // The number the button's own blurb prints: 11 s, against nominal's 1.8.
    ASSERT_TRUE(s.settle_time > 10.0 && s.settle_time < 13.0);
    // It does not wander further out than nominal on the way — sluggish, not
    // wild.  Both start 72.1 mm out and neither adds to it.
    ASSERT_TRUE(s.peak_radius < base.peak_radius + 0.001);
}

// "Aggressive visibly overshoots OR saturates."  It saturates, and it is worth
// being exact about why the first disjunct is not the one satisfied: on this
// plant the travel limits bite before the ball ever gets past the setpoint, so
// saturation LIMITS the approach rather than producing an overshoot.  Measured
// position overshoot past centre is under a millimetre at every tuning that
// keeps the ball, which is why the surface says "saturates" and not "overshoots".
void test_the_aggressive_preset_is_fast_and_saturates() {
    const SimResult s = runPreset(presetNamed("Aggressive"));
    const SimResult base = nominalRun();

    ASSERT_TRUE(!s.left_plate);
    ASSERT_TRUE(s.saturated);                    // and nominal, above, does not
    ASSERT_TRUE(s.settle_time < 0.5 * base.settle_time);
    ASSERT_TRUE(s.final_radius < 0.005);
    // The number the button's own blurb prints: 0.6 s.
    ASSERT_TRUE(s.settle_time > 0.4 && s.settle_time < 0.8);
    // CONTEXT.md points a visitor at BOTH badges in Plate Control, so both
    // have to actually light up.  They are different facts: the legs stop at
    // their travel limits, and separately the triple they were asked for has
    // no assembly at all.
    ASSERT_TRUE(s.clipped);
}

// Why the Aggressive preset is Q = 150 and not the 300 the ticket's own
// comment proposed.  That comment predates #23, which let the ball leave the
// plate: a tuning that slams the legs moves the plate out from under the ball
// rather than tilting under it, and at 300 a nudge throws the ball clean off.
//
// The honest statement of where the line is drawn is a comparison, not an
// absolute — Aggressive is offered because it is no more fragile than the
// tuning the demo already ships, against the hardest shove the interface can
// deliver.  Both DO lose the ball above `kMaxNudgeSpeed`; that is a property
// of the plate, not of the preset, and is why the slider stops where it does.
//
// Pinned here rather than left as prose because prose does not fail when
// somebody raises the slider's ceiling.  Twenty-four directions at the top of
// the slider: a shove the loop cannot reject is one the demo must not offer,
// wherever it points.
void test_aggressive_is_no_more_fragile_than_the_shipped_tuning() {
    const auto models = getBuiltinModels();
    const auto& e = cascadeModel(models);
    const Eigen::MatrixXd agg = gainForPreset(e, presetNamed("Aggressive"));
    const Eigen::MatrixXd nom = gainForPreset(e, nominalPreset());

    for (int i = 0; i < 24; ++i) {
        const double th = i * M_PI / 12.0;
        const Eigen::Vector4d kick(0.0, 0.0,
                                   kMaxNudgeSpeed * std::cos(th),
                                   kMaxNudgeSpeed * std::sin(th));
        ASSERT_TRUE(!runClosedLoop(e, nom, kick, 0.0, 0.0, 8.0).left_plate);
        ASSERT_TRUE(!runClosedLoop(e, agg, kick, 0.0, 0.0, 8.0).left_plate);
    }
}

// Whatever else a preset does, it has to leave a working demo behind: a
// visitor who clicks one and walks away must not come back to a ball on the
// floor.  Cheap to state, and the one property all three share.
void test_every_preset_brings_the_ball_home() {
    for (const auto& p : lqrPresets()) {
        const SimResult s = runPreset(p);
        ASSERT_TRUE(!s.left_plate);
        ASSERT_TRUE(s.settle_time > 0.0);
        ASSERT_TRUE(s.final_radius < 0.005);
    }
}

}  // namespace

int main() {
    test_state_vector_layout();
    test_state_deviation_follows_the_design_home();
    test_shape_gate();
    test_plant_identity_survives_the_float_round_trip();
    test_zero_error_commands_the_home_pose();
    test_command_is_home_minus_k_error();
    test_setpoint_moves_the_target();
    test_the_reference_velocity_enters_as_a_velocity_error();
    test_the_ball_velocity_weights_set_the_feedforwards_strength();
    test_command_clamps_to_servo_travel_then_to_the_workspace();
    test_a_reachable_clamp_is_left_on_the_stop();
    test_wrong_shape_commands_the_home_pose();
    test_servo_lag_is_first_order();
    test_default_weights_balance_the_ball();
    test_default_weights_reject_a_nudge();
    test_setpoint_is_tracked();
    test_looser_weights_are_visibly_sluggish();
    test_poor_tuning_visibly_degrades();
    test_unit_weights_park_in_the_friction_dead_band();
    test_the_nominal_preset_is_the_startup_tuning();
    test_the_presets_differ_only_in_the_ball_weights();
    test_the_lit_button_follows_the_weights();
    test_the_nominal_preset_settles_without_touching_the_stops();
    test_the_detuned_preset_is_sluggish_but_never_unstable();
    test_the_aggressive_preset_is_fast_and_saturates();
    test_aggressive_is_no_more_fragile_than_the_shipped_tuning();
    test_every_preset_brings_the_ball_home();
    std::printf("test_auto_balance: all passed\n");
    return 0;
}
