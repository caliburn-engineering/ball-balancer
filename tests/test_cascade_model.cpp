// tests/test_cascade_model.cpp
//
// The cascaded Ball-Balancer plant spans three subsystems that were previously
// modelled apart: servo lag, the 3-RRS mechanism's kinematics, and the rolling
// ball.  Two things can go wrong silently.  The mechanism block can drift from
// the kinematics the 3D view actually draws, and the ball block can pick up the
// wrong tilt sign — the same confusion `test_ball_sim` pins for the simulator,
// now in the linear model the controller is designed against.
#include "analysis/model_library.h"
#include "analysis/lqr.h"
#include "analysis/system_properties.h"
#include "table_kinematics.h"
#include "test_helpers.h"

#include <cmath>
#include <string>

using namespace caliburn;

namespace {

constexpr double kDeg = M_PI / 180.0;
constexpr double kRollingFactor = 5.0 / 7.0;

const ModelEntry& cascadeModel(const std::vector<ModelEntry>& models) {
    for (const auto& m : models)
        if (m.name.rfind("Ball-Balancer Cascade", 0) == 0) return m;
    std::fprintf(stderr, "FAIL: no cascade model in the library\n");
    std::exit(1);
}

double paramValue(const ModelEntry& e, const std::string& symbol) {
    for (const auto& p : e.params)
        if (p.symbol == symbol) return p.value;
    std::fprintf(stderr, "FAIL: no parameter '%s'\n", symbol.c_str());
    std::exit(1);
}

// The kinematics the model claims to be linearised about, rebuilt from the
// model's own published parameters — so this cannot drift from the builder.
Eigen::Matrix3d jacobianFor(const ModelEntry& e) {
    TableParams tp;
    tp.R_ground  = paramValue(e, "Rg");
    tp.R_table   = paramValue(e, "Rt");
    tp.L1        = paramValue(e, "L1");
    tp.L2        = paramValue(e, "L2");
    tp.alpha_min = 10.0 * kDeg;
    tp.alpha_max = 80.0 * kDeg;
    TableKinematics tk(tp);

    const double a_home = paramValue(e, "a0") * kDeg;
    const TablePose home = tk.home_pose(a_home);
    return tk.velocity_jacobian({a_home, a_home, a_home}, home);
}

// --- Shape: 3 leg commands in, roll/pitch/heave + ball x/y out ---
void test_dimensions() {
    const auto models = getBuiltinModels();
    const auto& e = cascadeModel(models);

    ASSERT_EQ(e.system.states(), 7);   // 3 servo + 4 ball
    ASSERT_EQ(e.system.inputs(), 3);   // three leg angles
    ASSERT_EQ(e.system.outputs(), 5);  // roll, pitch, heave, x, y
}

// --- The servo block is a first-order lag per leg, and nothing else ---
//
// "Nothing else" now carries a second job.  The plant's servos gained a rate
// limit in #32, and the limit deliberately did NOT come with it: it is a
// saturation outside the design model, accepted and documented exactly as the
// travel clamp is, not a redesign trigger.  `A(i,i) = -1/tau` with nothing
// beside it is what says the linearised cascade model — and so K, which is
// designed against it — is unchanged.
void test_servo_block() {
    const auto models = getBuiltinModels();
    const auto& e = cascadeModel(models);
    const double tau = paramValue(e, "\xcf\x84");  // tau

    for (int i = 0; i < 3; ++i) {
        ASSERT_NEAR(e.system.A(i, i), -1.0 / tau, 1e-12);
        ASSERT_NEAR(e.system.B(i, i), 1.0 / tau, 1e-12);
        for (int j = 0; j < 3; ++j)
            if (i != j) {
                ASSERT_NEAR(e.system.A(i, j), 0.0, 1e-12);
                ASSERT_NEAR(e.system.B(i, j), 0.0, 1e-12);
            }
        // A leg command drives its own servo and nothing else directly.
        for (int r = 3; r < 7; ++r)
            ASSERT_NEAR(e.system.B(r, i), 0.0, 1e-12);
    }
}

// --- The mechanism block IS the kinematics, not a copy of it ---
void test_pose_outputs_are_the_jacobian() {
    const auto models = getBuiltinModels();
    const auto& e = cascadeModel(models);
    const Eigen::Matrix3d Jv = jacobianFor(e);

    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            ASSERT_NEAR(e.system.C(i, j), Jv(i, j), 1e-9);

    // Ball position is read straight off the ball states.
    ASSERT_NEAR(e.system.C(3, 3), 1.0, 1e-12);
    ASSERT_NEAR(e.system.C(4, 4), 1.0, 1e-12);
    ASSERT_NEAR(e.system.D.cwiseAbs().maxCoeff(), 0.0, 1e-12);
}

// --- A pure pitch drives the ball along +x and leaves y alone ---
void test_pitch_drives_ball_in_positive_x() {
    const auto models = getBuiltinModels();
    const auto& e = cascadeModel(models);
    const double g = paramValue(e, "g");
    const Eigen::Matrix3d Jv = jacobianFor(e);

    // The leg deflection that produces one radian of pitch and no roll/heave.
    const Eigen::Vector3d d_alpha = Jv.colPivHouseholderQr().solve(
        Eigen::Vector3d(0.0, 1.0, 0.0));

    const Eigen::Vector3d a_vx = e.system.A.block(5, 0, 1, 3).transpose();
    const Eigen::Vector3d a_vy = e.system.A.block(6, 0, 1, 3).transpose();

    ASSERT_NEAR(a_vx.dot(d_alpha), kRollingFactor * g, 1e-6);
    ASSERT_NEAR(a_vy.dot(d_alpha), 0.0, 1e-6);
}

// --- A pure roll drives the ball along -y and leaves x alone ---
void test_roll_drives_ball_in_negative_y() {
    const auto models = getBuiltinModels();
    const auto& e = cascadeModel(models);
    const double g = paramValue(e, "g");
    const Eigen::Matrix3d Jv = jacobianFor(e);

    const Eigen::Vector3d d_alpha = Jv.colPivHouseholderQr().solve(
        Eigen::Vector3d(1.0, 0.0, 0.0));

    const Eigen::Vector3d a_vx = e.system.A.block(5, 0, 1, 3).transpose();
    const Eigen::Vector3d a_vy = e.system.A.block(6, 0, 1, 3).transpose();

    ASSERT_NEAR(a_vy.dot(d_alpha), -kRollingFactor * g, 1e-6);
    ASSERT_NEAR(a_vx.dot(d_alpha), 0.0, 1e-6);
}

// --- The ball integrates its own velocity, and heave never reaches it ---
void test_ball_kinematic_chain() {
    const auto models = getBuiltinModels();
    const auto& e = cascadeModel(models);
    const Eigen::Matrix3d Jv = jacobianFor(e);

    ASSERT_NEAR(e.system.A(3, 5), 1.0, 1e-12);  // x' = vx
    ASSERT_NEAR(e.system.A(4, 6), 1.0, 1e-12);  // y' = vy

    // Pure heave: the plate rises level, so the ball must not accelerate.
    const Eigen::Vector3d d_alpha = Jv.colPivHouseholderQr().solve(
        Eigen::Vector3d(0.0, 0.0, 1.0));
    const Eigen::Vector3d a_vx = e.system.A.block(5, 0, 1, 3).transpose();
    const Eigen::Vector3d a_vy = e.system.A.block(6, 0, 1, 3).transpose();
    ASSERT_NEAR(a_vx.dot(d_alpha), 0.0, 1e-6);
    ASSERT_NEAR(a_vy.dot(d_alpha), 0.0, 1e-6);
}

// --- Three legs are enough to command all seven states ---
void test_controllable() {
    const auto models = getBuiltinModels();
    const auto& e = cascadeModel(models);
    const PropertyResult ctrb = checkControllability(e.system);
    ASSERT_TRUE(ctrb.pass);
    ASSERT_EQ(ctrb.rank, 7);
}

// --- LQR closes the loop: the whole point of the cascade ---
void test_lqr_stabilises_the_cascade() {
    const auto models = getBuiltinModels();
    const auto& e = cascadeModel(models);

    Eigen::MatrixXd Q = Eigen::MatrixXd::Identity(7, 7);
    Q(3, 3) = 100.0;  // ball position is what we care about
    Q(4, 4) = 100.0;
    const Eigen::MatrixXd R = Eigen::MatrixXd::Identity(3, 3);

    const LqrResult res = computeLQR(e.system, Q, R);
    ASSERT_TRUE(res.success);
    ASSERT_EQ((int)res.K.rows(), 3);
    ASSERT_EQ((int)res.K.cols(), 7);

    // Every closed-loop pole strictly in the left half plane — the open-loop
    // plant has four poles at the origin, so this is the substantive claim.
    for (const auto& p : res.closed_loop_poles)
        ASSERT_TRUE(p.real() < -1e-9);
}

}  // namespace

int main() {
    test_dimensions();
    test_servo_block();
    test_pose_outputs_are_the_jacobian();
    test_pitch_drives_ball_in_positive_x();
    test_roll_drives_ball_in_negative_y();
    test_ball_kinematic_chain();
    test_controllable();
    test_lqr_stabilises_the_cascade();
    std::printf("test_cascade_model: all passed\n");
    return 0;
}
