// tests/test_lqr.cpp
#include "analysis/lqr.h"
#include "analysis/model_library.h"
#include "analysis/system_connect.h"
#include "test_helpers.h"

#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <Eigen/LU>
#include <algorithm>
#include <complex>
#include <cmath>
#include <cstdio>
#include <vector>

using namespace caliburn;
using cd = std::complex<double>;

namespace {

// The design states the oracle as K ~ [0.414, 0.318] with closed-loop poles at
// -0.86 +/- 0.86j.  Those two numbers are inconsistent with each other and with
// the plant: K = [0.414, 0.318] places the poles at -0.859 +/- 0.822j, while
// poles at -0.86 +/- 0.86j require K = [0.479, 0.320].  Neither is the CARE
// solution.
//
// Solving the 2x2 Riccati equation by hand for A = [0 1; -1 -1.4], B = [0;1],
// Q = I, R = 1 gives P12 = sqrt(2) - 1 and P22 = -1.4 + sqrt(0.96 + 2*sqrt(2)),
// so K = [0.4142135624, 0.5463882256]; the stable eigenvalues of the
// Hamiltonian, computed independently of any Riccati solver, are
// -0.9731941128 +/- 0.6834521060j.  Three derivations agree, so the design's
// second gain component and its pole pair are simply wrong, and the values
// below are the oracle.  The residual assertion is what makes that call safe
// to make: it holds against the equation itself, not a transcribed number.
constexpr double kK1 = 0.4142135624;
constexpr double kK2 = 0.5463882256;
constexpr double kPoleRe = -0.9731941128;
constexpr double kPoleIm = 0.6834521060;

// The same plant with R = 4 instead of 1.  Every assertion above holds with
// R = I, where R^-1 is the identity and a solver that ignored R entirely would
// pass -- a mutation that drops R^-1 from K = R^-1B'P survives the whole suite
// without this case.  Hand-solved the same way: P12 = -R + sqrt(R^2 + R) and
// P22 = -1.4R + sqrt(1.96R^2 + R(1 + 2*P12)), with K = P_row / R.
constexpr double kRHeavy = 4.0;
constexpr double kRHeavyK1 = 0.1180339887;
constexpr double kRHeavyK2 = 0.1639910414;
constexpr double kRHeavyPoleRe = -0.7819955207;
constexpr double kRHeavyPoleIm = 0.7117000733;

constexpr double kTol = 1e-6;       // gains and poles
constexpr double kResTol = 1e-8;    // Riccati residual, per the design

Eigen::MatrixXd eye(int n) { return Eigen::MatrixXd::Identity(n, n); }

LinearSystem massSpringDamper() {
    LinearSystem sys;
    sys.A = Eigen::MatrixXd(2, 2);
    sys.A << 0.0, 1.0,
            -1.0, -1.4;
    sys.B = Eigen::MatrixXd(2, 1);
    sys.B << 0.0, 1.0;
    sys.C = eye(2);
    sys.D = Eigen::MatrixXd::Zero(2, 1);
    return sys;
}

// The Ball-Balancer plate: two decoupled double integrators, cross-coupled
// input to output (tilt about y drives x).  Two inputs, so it exercises the
// m > 1 path -- R is a matrix to invert, not a scalar to divide by.
LinearSystem ballPlate() {
    constexpr double kg = (5.0 / 7.0) * 9.81;
    LinearSystem sys;
    sys.A = Eigen::MatrixXd::Zero(4, 4);
    sys.A(0, 2) = 1.0;
    sys.A(1, 3) = 1.0;
    sys.B = Eigen::MatrixXd::Zero(4, 2);
    sys.B(2, 1) = kg;
    sys.B(3, 0) = kg;
    sys.C = Eigen::MatrixXd::Zero(2, 4);
    sys.C(0, 0) = 1.0;
    sys.C(1, 1) = 1.0;
    sys.D = Eigen::MatrixXd::Zero(2, 2);
    return sys;
}

// The assertion that earns its place: a solver can return plausible gains while
// satisfying no equation at all, and only this catches that.
double riccatiResidual(const LinearSystem& sys, const Eigen::MatrixXd& Q,
                       const Eigen::MatrixXd& R, const Eigen::MatrixXd& P) {
    const Eigen::MatrixXd& A = sys.A;
    const Eigen::MatrixXd& B = sys.B;
    const Eigen::MatrixXd res = A.transpose() * P + P * A
                              - P * B * R.inverse() * B.transpose() * P + Q;
    return res.norm();
}

std::vector<cd> sortedPoles(std::vector<cd> poles) {
    std::sort(poles.begin(), poles.end(), [](const cd& a, const cd& b) {
        if (a.real() != b.real()) return a.real() < b.real();
        return a.imag() < b.imag();
    });
    return poles;
}

// --- The oracle -------------------------------------------------------------

void test_mass_spring_damper_oracle() {
    const LinearSystem sys = massSpringDamper();
    const Eigen::MatrixXd Q = eye(2);
    const Eigen::MatrixXd R = eye(1);

    const LqrResult r = computeLQR(sys, Q, R);
    ASSERT_TRUE(r.success);
    ASSERT_TRUE(r.error.empty());

    ASSERT_EQ(r.K.rows(), 1);
    ASSERT_EQ(r.K.cols(), 2);
    ASSERT_NEAR(r.K(0, 0), kK1, kTol);
    ASSERT_NEAR(r.K(0, 1), kK2, kTol);

    const std::vector<cd> poles = sortedPoles(r.closed_loop_poles);
    ASSERT_EQ(poles.size(), 2u);
    ASSERT_NEAR(poles[0].real(), kPoleRe, kTol);
    ASSERT_NEAR(poles[0].imag(), -kPoleIm, kTol);
    ASSERT_NEAR(poles[1].real(), kPoleRe, kTol);
    ASSERT_NEAR(poles[1].imag(), kPoleIm, kTol);

    ASSERT_TRUE(riccatiResidual(sys, Q, R, r.P) < kResTol);
}

void test_riccati_solution_is_symmetric_positive_definite() {
    const LinearSystem sys = massSpringDamper();
    const LqrResult r = computeLQR(sys, eye(2),
                                   eye(1));
    ASSERT_TRUE(r.success);
    ASSERT_TRUE((r.P - r.P.transpose()).norm() < 1e-12);
    ASSERT_TRUE(Eigen::LLT<Eigen::MatrixXd>(r.P).info() == Eigen::Success);
}

// K = R^-1 B'P is the definition, and it is cheap to check that the returned
// pair is internally consistent rather than two independently wrong matrices.
void test_gain_matches_riccati_solution() {
    const LinearSystem sys = massSpringDamper();
    const Eigen::MatrixXd Q = eye(2);
    const Eigen::MatrixXd R = eye(1);
    const LqrResult r = computeLQR(sys, Q, R);
    ASSERT_TRUE(r.success);
    const Eigen::MatrixXd expected = R.inverse() * sys.B.transpose() * r.P;
    ASSERT_TRUE((r.K - expected).norm() < 1e-12);
}

// The poles the solver reports must be the poles of the system it actually
// closes -- the existing stateFeedbackClose is the arbiter, not a second
// eigenvalue computation inside the solver.
void test_closed_loop_poles_match_state_feedback_close() {
    const LinearSystem sys = massSpringDamper();
    const LqrResult r = computeLQR(sys, eye(2),
                                   eye(1));
    ASSERT_TRUE(r.success);

    const LinearSystem cl = stateFeedbackClose(sys, r.K);
    Eigen::EigenSolver<Eigen::MatrixXd> es(cl.A);
    std::vector<cd> actual(es.eigenvalues().data(),
                           es.eigenvalues().data() + es.eigenvalues().size());
    actual = sortedPoles(actual);
    const std::vector<cd> reported = sortedPoles(r.closed_loop_poles);
    ASSERT_EQ(actual.size(), reported.size());
    for (size_t i = 0; i < actual.size(); ++i) {
        ASSERT_NEAR(actual[i].real(), reported[i].real(), 1e-12);
        ASSERT_NEAR(actual[i].imag(), reported[i].imag(), 1e-12);
    }
}

// R is the only weight that enters the gain through an inverse, so it is the
// one a solver can silently drop.  This is the same oracle as above, derived
// the same way, at a weight where R^-1 is not the identity.
void test_input_weight_scales_the_gain() {
    const LinearSystem sys = massSpringDamper();
    const Eigen::MatrixXd Q = eye(2);
    Eigen::MatrixXd R(1, 1);
    R << kRHeavy;

    const LqrResult r = computeLQR(sys, Q, R);
    ASSERT_TRUE(r.success);
    ASSERT_NEAR(r.K(0, 0), kRHeavyK1, kTol);
    ASSERT_NEAR(r.K(0, 1), kRHeavyK2, kTol);

    const std::vector<cd> poles = sortedPoles(r.closed_loop_poles);
    ASSERT_NEAR(poles[0].real(), kRHeavyPoleRe, kTol);
    ASSERT_NEAR(poles[0].imag(), -kRHeavyPoleIm, kTol);
    ASSERT_NEAR(poles[1].imag(), kRHeavyPoleIm, kTol);

    ASSERT_TRUE(riccatiResidual(sys, Q, R, r.P) < kResTol);

    // K = R^-1B'P must hold at R != I too, not just where the inverse vanishes.
    const Eigen::MatrixXd expected = R.inverse() * sys.B.transpose() * r.P;
    ASSERT_TRUE((r.K - expected).norm() < 1e-12);
}

// --- Weighting behaviour ----------------------------------------------------

// Heavier state weighting buys tighter regulation with more effort: the
// closed-loop poles move left.  This is the property the UI sliders exist to
// demonstrate, and it is direction, not magnitude, so no oracle is needed.
void test_heavier_q_moves_poles_left() {
    const LinearSystem sys = massSpringDamper();
    const Eigen::MatrixXd R = eye(1);

    const LqrResult light = computeLQR(sys, eye(2), R);
    const LqrResult heavy =
        computeLQR(sys, 100.0 * eye(2), R);
    ASSERT_TRUE(light.success);
    ASSERT_TRUE(heavy.success);

    double light_max = -1e300, heavy_max = -1e300;
    for (const cd& p : light.closed_loop_poles)
        light_max = std::max(light_max, p.real());
    for (const cd& p : heavy.closed_loop_poles)
        heavy_max = std::max(heavy_max, p.real());
    ASSERT_TRUE(heavy_max < light_max);
}

// Penalising effort buys a quieter actuator and a slower loop: the gain shrinks
// and the poles move right.  The mirror of the Q test, and the other half of
// what the design's Q/R sliders exist to show.
void test_heavier_r_shrinks_the_gain() {
    const LinearSystem sys = massSpringDamper();
    const Eigen::MatrixXd Q = eye(2);

    Eigen::MatrixXd R_light(1, 1), R_heavy(1, 1);
    R_light << 1.0;
    R_heavy << 100.0;

    const LqrResult light = computeLQR(sys, Q, R_light);
    const LqrResult heavy = computeLQR(sys, Q, R_heavy);
    ASSERT_TRUE(light.success);
    ASSERT_TRUE(heavy.success);

    ASSERT_TRUE(heavy.K.norm() < light.K.norm());

    double light_max = -1e300, heavy_max = -1e300;
    for (const cd& p : light.closed_loop_poles)
        light_max = std::max(light_max, p.real());
    for (const cd& p : heavy.closed_loop_poles)
        heavy_max = std::max(heavy_max, p.real());
    ASSERT_TRUE(heavy_max > light_max);
}

// --- MIMO -------------------------------------------------------------------

void test_ball_plate_two_input() {
    const LinearSystem sys = ballPlate();
    Eigen::MatrixXd Q = eye(4);
    Q(0, 0) = 50.0;
    Q(1, 1) = 50.0;
    const Eigen::MatrixXd R = eye(2);

    const LqrResult r = computeLQR(sys, Q, R);
    ASSERT_TRUE(r.success);
    ASSERT_EQ(r.K.rows(), 2);
    ASSERT_EQ(r.K.cols(), 4);

    ASSERT_EQ(r.closed_loop_poles.size(), 4u);
    for (const cd& p : r.closed_loop_poles) ASSERT_TRUE(p.real() < 0.0);

    ASSERT_TRUE(riccatiResidual(sys, Q, R, r.P) < kResTol);
}

// The plate's two axes are independent and identically weighted, so the design
// must come out symmetric under swapping them.  A solver that special-cases the
// first input, or that leaks a column index, breaks this.
void test_ball_plate_axes_are_symmetric() {
    const LinearSystem sys = ballPlate();
    const LqrResult r = computeLQR(sys, eye(4),
                                   eye(2));
    ASSERT_TRUE(r.success);
    // input 0 drives state 3 (y-chain), input 1 drives state 2 (x-chain)
    ASSERT_NEAR(r.K(0, 1), r.K(1, 0), 1e-9);  // position feedback
    ASSERT_NEAR(r.K(0, 3), r.K(1, 2), 1e-9);  // velocity feedback
}

// Two inputs weighted differently: the expensive axis must end up with the
// smaller gain.  This is the m > 1 path with a genuine matrix inverse in it,
// where R is neither a scalar nor the identity.
void test_ball_plate_with_unequal_input_weights() {
    const LinearSystem sys = ballPlate();
    const Eigen::MatrixXd Q = eye(4);
    Eigen::MatrixXd R = eye(2);
    R(1, 1) = 100.0;  // input 1 (drives the x-chain) is expensive

    const LqrResult r = computeLQR(sys, Q, R);
    ASSERT_TRUE(r.success);
    ASSERT_TRUE(riccatiResidual(sys, Q, R, r.P) < kResTol);

    // Input 1 drives state 2, input 0 drives state 3 (see ballPlate).  With
    // input 1 penalised 100x, its position feedback must be the weaker of the
    // two -- the symmetry asserted above is broken in the direction R dictates.
    ASSERT_TRUE(std::abs(r.K(1, 0)) < std::abs(r.K(0, 1)));

    for (const cd& p : r.closed_loop_poles) ASSERT_TRUE(p.real() < 0.0);
}

// --- Refinement and residual ------------------------------------------------

// Every passing case must achieve a relative Riccati residual below 1e-12
// after Newton/Kleinman refinement.  This is tested per-system so a future
// regression is traceable to the plant that regressed.
void test_residual_below_machine_precision_mass_spring_damper() {
    const LinearSystem sys = massSpringDamper();
    const LqrResult r = computeLQR(sys, eye(2), eye(1));
    ASSERT_TRUE(r.success);
    ASSERT_TRUE(r.residual < 1e-12);
}

void test_residual_below_machine_precision_r_heavy() {
    const LinearSystem sys = massSpringDamper();
    Eigen::MatrixXd R(1, 1);
    R << kRHeavy;
    const LqrResult r = computeLQR(sys, eye(2), R);
    ASSERT_TRUE(r.success);
    ASSERT_TRUE(r.residual < 1e-12);
}

void test_residual_below_machine_precision_ball_plate() {
    const LinearSystem sys = ballPlate();
    Eigen::MatrixXd Q = eye(4);
    Q(0, 0) = 50.0;
    Q(1, 1) = 50.0;
    const LqrResult r = computeLQR(sys, Q, eye(2));
    ASSERT_TRUE(r.success);
    ASSERT_TRUE(r.residual < 1e-12);
}

// Reports the pre-refinement relative residual on the cascade plant (7 states,
// 14×14 Hamiltonian) so the evidence for or against Hamiltonian balancing is
// on the record.  The assertion is on the post-refinement residual; the
// pre-refinement number is printed for that record.
void test_pre_refinement_residual_on_cascade_plant() {
    const std::vector<ModelEntry>& models = getBuiltinModels();
    const ModelEntry* cascade = nullptr;
    for (const auto& m : models)
        if (isCascadeModel(m)) { cascade = &m; break; }
    ASSERT_TRUE(cascade != nullptr);

    // Nominal preset weights: leg states = 1, ball position = 100, velocity = 10
    const int n = cascade->system.states();
    const int m = cascade->system.inputs();
    ASSERT_EQ(n, 7);
    ASSERT_EQ(m, 3);
    Eigen::VectorXd q_diag = Eigen::VectorXd::Ones(n);
    q_diag(3) = q_diag(4) = 100.0;  // ball position
    q_diag(5) = q_diag(6) = 10.0;   // ball velocity
    const Eigen::MatrixXd Q = q_diag.asDiagonal().toDenseMatrix();
    const Eigen::MatrixXd R = Eigen::MatrixXd::Identity(m, m);

    const LqrResult r = computeLQR(cascade->system, Q, R);
    ASSERT_TRUE(r.success);

    std::printf("  cascade plant pre-refinement relative residual:  %.3e\n",
                r.pre_refinement_residual);
    std::printf("  cascade plant post-refinement relative residual: %.3e\n",
                r.residual);

    ASSERT_TRUE(r.residual < 1e-12);
}

// --- Rejected inputs --------------------------------------------------------

void test_rejects_non_positive_definite_r() {
    const LinearSystem sys = massSpringDamper();
    Eigen::MatrixXd R(1, 1);
    R << 0.0;
    const LqrResult r = computeLQR(sys, eye(2), R);
    ASSERT_TRUE(!r.success);
    ASSERT_TRUE(!r.error.empty());
    ASSERT_EQ(r.K.size(), 0);
}

void test_rejects_negative_r() {
    const LinearSystem sys = massSpringDamper();
    Eigen::MatrixXd R(1, 1);
    R << -1.0;
    const LqrResult r = computeLQR(sys, eye(2), R);
    ASSERT_TRUE(!r.success);
}

void test_rejects_wrong_sized_q() {
    const LinearSystem sys = massSpringDamper();
    const LqrResult r = computeLQR(sys, eye(3),
                                   eye(1));
    ASSERT_TRUE(!r.success);
    ASSERT_TRUE(!r.error.empty());
}

void test_rejects_wrong_sized_r() {
    const LinearSystem sys = massSpringDamper();
    const LqrResult r = computeLQR(sys, eye(2),
                                   eye(2));
    ASSERT_TRUE(!r.success);
    ASSERT_TRUE(!r.error.empty());
}

void test_rejects_asymmetric_q() {
    const LinearSystem sys = massSpringDamper();
    Eigen::MatrixXd Q(2, 2);
    Q << 1.0, 2.0,
         0.0, 1.0;
    const LqrResult r = computeLQR(sys, Q, eye(1));
    ASSERT_TRUE(!r.success);
}

void test_rejects_indefinite_q() {
    const LinearSystem sys = massSpringDamper();
    Eigen::MatrixXd Q(2, 2);
    Q << -1.0, 0.0,
          0.0, 1.0;
    const LqrResult r = computeLQR(sys, Q, eye(1));
    ASSERT_TRUE(!r.success);
}

// A stabilizable-but-not-controllable plant: the uncontrollable mode is
// stable at -2, so it decays on its own.  The CARE has a solution and the
// solver must accept it.
void test_accepts_stabilizable_but_uncontrollable_plant() {
    // state 0: x0' = -x0 + u  (controllable, stable at -1)
    // state 1: x1' = -2*x1    (uncontrollable, stable at -2)
    LinearSystem sys;
    sys.A = Eigen::MatrixXd::Zero(2, 2);
    sys.A(0, 0) = -1.0;
    sys.A(1, 1) = -2.0;
    sys.B = Eigen::MatrixXd::Zero(2, 1);
    sys.B(0, 0) = 1.0;
    sys.C = eye(2);
    sys.D = Eigen::MatrixXd::Zero(2, 1);

    // Zero weight on the uncontrollable state keeps the cost integral finite.
    Eigen::MatrixXd Q = eye(2);
    Q(1, 1) = 0.0;

    const LqrResult r = computeLQR(sys, Q, eye(1));
    ASSERT_TRUE(r.success);
    ASSERT_TRUE(r.error.empty());
    ASSERT_EQ(r.K.rows(), 1);
    ASSERT_EQ(r.K.cols(), 2);
    for (const cd& p : r.closed_loop_poles) ASSERT_TRUE(p.real() < 0.0);
}

// An uncontrollable, unstabilizable plant has no finite-cost solution.  The
// second state is disconnected from the input and unstable, so no gain can
// move it.
void test_rejects_uncontrollable_plant() {
    LinearSystem sys;
    sys.A = Eigen::MatrixXd::Zero(2, 2);
    sys.A(0, 0) = -1.0;
    sys.A(1, 1) = 1.0;
    sys.B = Eigen::MatrixXd::Zero(2, 1);
    sys.B(0, 0) = 1.0;
    sys.C = eye(2);
    sys.D = Eigen::MatrixXd::Zero(2, 1);

    const LqrResult r = computeLQR(sys, eye(2),
                                   eye(1));
    ASSERT_TRUE(!r.success);
    ASSERT_TRUE(!r.error.empty());
}

void test_rejects_empty_system() {
    LinearSystem sys;
    const LqrResult r = computeLQR(sys, Eigen::MatrixXd(), Eigen::MatrixXd());
    ASSERT_TRUE(!r.success);
}

// --- Residual gate ---------------------------------------------------------
//
// The gate lives at 1e-8.  Normal solutions reach 1e-12 to 1e-15 after
// Kleinman refinement (pinned above), so the gate only trips on genuine
// numerical breakage — a sign-function result that converged to the wrong
// invariant subspace, or a P extraction that was numerically poor.  Such
// cases are not reachable from standard well-conditioned inputs; the
// following tests verify the gap between the suite's 1e-12 assertion and the
// 1e-8 gate (AC3), and that the gate's error message is distinguishable from
// input-rejection failures (AC2).

// AC3: every passing case in the suite lands well below the gate.
// (Each residual test above already asserts < 1e-12; this one makes the gap
// explicit and is the single place that would fail first if a regression
// pushed any solve above 1e-8.)
void test_residual_gate_is_never_tripped_by_suite_inputs() {
    struct Case { const char* name; LinearSystem sys; Eigen::MatrixXd Q, R; };

    const LinearSystem msd = massSpringDamper();
    const LinearSystem bp  = ballPlate();

    Eigen::MatrixXd R_heavy(1, 1);
    R_heavy << kRHeavy;

    Eigen::MatrixXd Q_bp = eye(4);
    Q_bp(0, 0) = 50.0; Q_bp(1, 1) = 50.0;

    const Case cases[] = {
        {"mass_spring_damper",      msd, eye(2), eye(1)},
        {"mass_spring_damper_Rhvy", msd, eye(2), R_heavy},
        {"ball_plate",              bp,  Q_bp,   eye(2)},
    };
    for (const auto& c : cases) {
        const LqrResult r = computeLQR(c.sys, c.Q, c.R);
        ASSERT_TRUE(r.success);
        // Well below the 1e-8 gate: genuine breakage would not look like this.
        ASSERT_TRUE(r.residual < kLqrResidualGate);
    }
}

// AC2: input-rejection and solver-convergence failures do not say "residual",
// establishing that a gate failure (which does) is distinguishable from them.
void test_input_rejection_messages_do_not_contain_residual() {
    using R1 = LqrResult;
    const LinearSystem sys = massSpringDamper();

    // Wrong-size Q
    const R1 r_qsize = computeLQR(sys, eye(3), eye(1));
    ASSERT_TRUE(!r_qsize.success);
    ASSERT_TRUE(r_qsize.error.find("residual") == std::string::npos);

    // Wrong-size R
    const R1 r_rsize = computeLQR(sys, eye(2), eye(2));
    ASSERT_TRUE(!r_rsize.success);
    ASSERT_TRUE(r_rsize.error.find("residual") == std::string::npos);

    // Asymmetric Q
    Eigen::MatrixXd Qasym(2, 2);
    Qasym << 1.0, 2.0, 0.0, 1.0;
    const R1 r_qasym = computeLQR(sys, Qasym, eye(1));
    ASSERT_TRUE(!r_qasym.success);
    ASSERT_TRUE(r_qasym.error.find("residual") == std::string::npos);

    // Non-positive-definite R
    const R1 r_rpd = computeLQR(sys, eye(2), Eigen::MatrixXd::Zero(1, 1));
    ASSERT_TRUE(!r_rpd.success);
    ASSERT_TRUE(r_rpd.error.find("residual") == std::string::npos);

    // Not stabilizable
    LinearSystem unstab;
    unstab.A = Eigen::MatrixXd::Zero(2, 2);
    unstab.A(0, 0) = -1.0; unstab.A(1, 1) = 1.0;
    unstab.B = Eigen::MatrixXd::Zero(2, 1);
    unstab.B(0, 0) = 1.0;
    unstab.C = eye(2); unstab.D = Eigen::MatrixXd::Zero(2, 1);
    const R1 r_unstab = computeLQR(unstab, eye(2), eye(1));
    ASSERT_TRUE(!r_unstab.success);
    ASSERT_TRUE(r_unstab.error.find("residual") == std::string::npos);
}

// AC1: when the gate trips, the result is a failure with no gain, no Riccati
// solution and no poles, and the message states the residual.  No
// well-conditioned plant reaches the 1e-8 gate, so the test lowers it to zero:
// any solve with a nonzero residual then trips it through the same code path.
void test_residual_gate_failure_clears_result_and_states_residual() {
    const LinearSystem sys = ballPlate();
    Eigen::MatrixXd Q = eye(4);
    Q(0, 0) = 50.0;
    Q(1, 1) = 50.0;

    const LqrResult ok = computeLQR(sys, Q, eye(2));
    ASSERT_TRUE(ok.success);
    ASSERT_TRUE(ok.residual > 0.0);  // otherwise a zero gate could not trip

    const LqrResult r = computeLQR(sys, Q, eye(2), 0.0);
    ASSERT_TRUE(!r.success);
    ASSERT_EQ(r.K.size(), 0);
    ASSERT_EQ(r.P.size(), 0);
    ASSERT_TRUE(r.closed_loop_poles.empty());
    ASSERT_TRUE(r.error.find("Riccati residual") != std::string::npos);

    char residual_text[32];
    std::snprintf(residual_text, sizeof(residual_text), "%.3e", ok.residual);
    ASSERT_TRUE(r.error.find(residual_text) != std::string::npos);
}

}  // namespace

int main() {
    test_mass_spring_damper_oracle();
    test_riccati_solution_is_symmetric_positive_definite();
    test_gain_matches_riccati_solution();
    test_closed_loop_poles_match_state_feedback_close();
    test_input_weight_scales_the_gain();
    test_heavier_q_moves_poles_left();
    test_heavier_r_shrinks_the_gain();
    test_ball_plate_two_input();
    test_ball_plate_axes_are_symmetric();
    test_ball_plate_with_unequal_input_weights();
    test_residual_below_machine_precision_mass_spring_damper();
    test_residual_below_machine_precision_r_heavy();
    test_residual_below_machine_precision_ball_plate();
    test_pre_refinement_residual_on_cascade_plant();
    test_accepts_stabilizable_but_uncontrollable_plant();
    test_rejects_non_positive_definite_r();
    test_rejects_negative_r();
    test_rejects_wrong_sized_q();
    test_rejects_wrong_sized_r();
    test_rejects_asymmetric_q();
    test_rejects_indefinite_q();
    test_rejects_uncontrollable_plant();
    test_rejects_empty_system();
    test_residual_gate_is_never_tripped_by_suite_inputs();
    test_input_rejection_messages_do_not_contain_residual();
    test_residual_gate_failure_clears_result_and_states_residual();
    std::printf("All lqr tests passed.\n");
    return 0;
}
