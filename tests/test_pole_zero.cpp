// tests/test_pole_zero.cpp
#include "test_helpers.h"
#include "analysis/pole_zero.h"
#include <algorithm>
#include <cmath>
#include <cstdio>

// Helper: find the pole closest to a given complex value
std::complex<double> find_nearest(
    const std::vector<std::complex<double>>& vec,
    std::complex<double> target) {
    double best = 1e30;
    std::complex<double> result = vec[0];
    for (const auto& v : vec) {
        double d = std::abs(v - target);
        if (d < best) { best = d; result = v; }
    }
    return result;
}

void test_first_order_pole() {
    caliburn::LinearSystem sys;
    sys.A = (Eigen::MatrixXd(1, 1) << -1).finished();
    sys.B = (Eigen::MatrixXd(1, 1) << 1).finished();
    sys.C = (Eigen::MatrixXd(1, 1) << 1).finished();
    sys.D = Eigen::MatrixXd::Zero(1, 1);

    auto result = caliburn::computePoleZero(sys, 0, 0);
    ASSERT_EQ((int)result.poles.size(), 1);
    ASSERT_NEAR(result.poles[0].real(), -1.0, 1e-10);
    ASSERT_NEAR(result.poles[0].imag(), 0.0, 1e-10);
    ASSERT_TRUE(result.is_stable);
    ASSERT_TRUE(result.zeros.empty());
}

void test_second_order_complex_poles() {
    caliburn::LinearSystem sys;
    sys.A = (Eigen::MatrixXd(2, 2) << 0, 1, -1, -1.4).finished();
    sys.B = (Eigen::MatrixXd(2, 1) << 0, 1).finished();
    sys.C = (Eigen::MatrixXd(1, 2) << 1, 0).finished();
    sys.D = Eigen::MatrixXd::Zero(1, 1);

    auto result = caliburn::computePoleZero(sys, 0, 0);
    ASSERT_EQ((int)result.poles.size(), 2);

    auto p1 = find_nearest(result.poles, {-0.7, 0.7141});
    ASSERT_NEAR(p1.real(), -0.7, 1e-4);
    ASSERT_NEAR(std::abs(p1.imag()), 0.7141, 1e-3);
    ASSERT_TRUE(result.is_stable);
}

void test_unstable_system() {
    caliburn::LinearSystem sys;
    sys.A = Eigen::MatrixXd::Zero(4, 4);
    sys.A(0, 1) = 1;
    sys.A(1, 2) = -0.981;
    sys.A(2, 3) = 1;
    sys.A(3, 2) = 21.582;
    sys.B = (Eigen::MatrixXd(4, 1) << 0, 1, 0, -2).finished();
    sys.C = Eigen::MatrixXd::Zero(2, 4);
    sys.C(0, 0) = 1; sys.C(1, 2) = 1;
    sys.D = Eigen::MatrixXd::Zero(2, 1);

    auto result = caliburn::computePoleZero(sys, 0, 0);
    ASSERT_TRUE(!result.is_stable);

    auto p = find_nearest(result.poles, {4.645, 0});
    ASSERT_NEAR(p.real(), 4.645, 0.01);
}

void test_transmission_zeros() {
    // G(s) = (s+3)/((s+1)(s+2)), zero at s = -3
    caliburn::LinearSystem sys;
    sys.A = (Eigen::MatrixXd(2, 2) << 0, 1, -2, -3).finished();
    sys.B = (Eigen::MatrixXd(2, 1) << 0, 1).finished();
    sys.C = (Eigen::MatrixXd(1, 2) << 3, 1).finished();
    sys.D = Eigen::MatrixXd::Zero(1, 1);

    auto result = caliburn::computePoleZero(sys, 0, 0);

    auto p1 = find_nearest(result.poles, {-1, 0});
    auto p2 = find_nearest(result.poles, {-2, 0});
    ASSERT_NEAR(p1.real(), -1.0, 1e-6);
    ASSERT_NEAR(p2.real(), -2.0, 1e-6);

    ASSERT_EQ((int)result.zeros.size(), 1);
    ASSERT_NEAR(result.zeros[0].real(), -3.0, 1e-6);
    ASSERT_NEAR(result.zeros[0].imag(), 0.0, 1e-6);
}

void test_zeros_with_nonzero_D() {
    // G(s) = (s+1)/(s+2), zero at s = -1
    // State-space: A=[-2], B=[1], C=[-1], D=[1]
    caliburn::LinearSystem sys;
    sys.A = (Eigen::MatrixXd(1, 1) << -2).finished();
    sys.B = (Eigen::MatrixXd(1, 1) << 1).finished();
    sys.C = (Eigen::MatrixXd(1, 1) << -1).finished();
    sys.D = (Eigen::MatrixXd(1, 1) << 1).finished();

    auto result = caliburn::computePoleZero(sys, 0, 0);
    ASSERT_EQ((int)result.zeros.size(), 1);
    ASSERT_NEAR(result.zeros[0].real(), -1.0, 1e-6);
}

void test_root_locus_endpoints() {
    caliburn::LinearSystem sys;
    sys.A = (Eigen::MatrixXd(1, 1) << -1).finished();
    sys.B = (Eigen::MatrixXd(1, 1) << 1).finished();
    sys.C = (Eigen::MatrixXd(1, 1) << 1).finished();
    sys.D = Eigen::MatrixXd::Zero(1, 1);

    auto locus = caliburn::computeRootLocus(sys, 0, 0, 0.0, 50.0, 100);
    ASSERT_EQ((int)locus.size(), 100);

    ASSERT_NEAR(locus.front().poles[0].real(), -1.0, 1e-6);
    ASSERT_NEAR(locus.back().poles[0].real(), -51.0, 1e-6);
}

void test_state_feedback_locus() {
    caliburn::LinearSystem sys;
    sys.A = (Eigen::MatrixXd(2, 2) << 0, 1, -1, -1).finished();
    sys.B = (Eigen::MatrixXd(2, 1) << 0, 1).finished();
    sys.C = (Eigen::MatrixXd(1, 2) << 1, 0).finished();
    sys.D = Eigen::MatrixXd::Zero(1, 1);

    Eigen::MatrixXd K(1, 2);
    K << 1, 0;

    auto locus = caliburn::computeStateFeedbackLocus(sys, K, 0.0, 3.0, 50);
    ASSERT_EQ((int)locus.size(), 50);

    auto p0 = find_nearest(locus.front().poles, {-0.5, 0.866});
    ASSERT_NEAR(p0.real(), -0.5, 1e-3);
    ASSERT_NEAR(std::abs(p0.imag()), 0.866, 0.01);

    double expected_imag = std::sqrt(15.0) / 2.0;
    auto p3 = find_nearest(locus.back().poles, {-0.5, expected_imag});
    ASSERT_NEAR(p3.real(), -0.5, 1e-3);
    ASSERT_NEAR(std::abs(p3.imag()), expected_imag, 0.02);
}

void test_match_poles_uses_optimal_assignment() {
    // Three-pole scenario where greedy nearest-neighbour gives a suboptimal
    // assignment.  prev = {0, 0.3, -100}, curr (as returned by the
    // eigenvalue solver) = {0.2, -10, -100.1}.
    //
    // Greedy processes prev[0]=0 first and picks curr[0]=0.2 (nearest, dist
    // 0.2), leaving curr[1]=-10 for prev[1]=0.3 (dist 10.3).  Greedy total:
    // 0.2 + 10.3 + 0.1 = 10.6.
    //
    // Optimal: prev[0] -> -10 (dist 10), prev[1] -> 0.2 (dist 0.1),
    //          prev[2] -> -100.1 (dist 0.1).  Total: 10.2.
    std::vector<std::complex<double>> prev = {{0.0, 0}, {0.3, 0}, {-100.0, 0}};
    std::vector<std::complex<double>> curr = {{0.2, 0}, {-10.0, 0}, {-100.1, 0}};
    caliburn::matchPoles(curr, prev);
    ASSERT_NEAR(curr[0].real(), -10.0,  1e-10);
    ASSERT_NEAR(curr[1].real(),   0.2,  1e-10);
    ASSERT_NEAR(curr[2].real(), -100.1, 1e-10);
}

void test_root_locus_sweep_branch_identity() {
    // G(s) = k / ((s+1)(s+3)): poles start at -1 and -3.  Observable
    // canonical form: A = [[0,1],[-3,-4]], B = [[0],[1]], C = [[1,0]].
    //
    // As k increases from 0 to 20 the two branches merge into a complex
    // conjugate pair on Re = -2.  With optimal branch matching each labelled
    // branch moves continuously; the test checks that no branch jumps by
    // more than 1 unit between consecutive gain steps (step size 0.1).
    caliburn::LinearSystem sys;
    sys.A = (Eigen::MatrixXd(2, 2) << 0, 1, -3, -4).finished();
    sys.B = (Eigen::MatrixXd(2, 1) << 0, 1).finished();
    sys.C = (Eigen::MatrixXd(1, 2) << 1, 0).finished();
    sys.D = Eigen::MatrixXd::Zero(1, 1);

    auto locus = caliburn::computeRootLocus(sys, 0, 0, 0.0, 20.0, 201);
    ASSERT_EQ((int)locus.size(), 201);

    // Verify continuity: no branch jumps more than 1 unit per step.
    for (std::size_t i = 1; i < locus.size(); ++i) {
        for (int b = 0; b < 2; ++b) {
            double jump = std::abs(locus[i].poles[b] - locus[i - 1].poles[b]);
            ASSERT_TRUE(jump < 1.0);
        }
    }

    // At high gain both poles sit on Re = -2.
    for (int b = 0; b < 2; ++b)
        ASSERT_NEAR(locus.back().poles[b].real(), -2.0, 0.1);
}

void test_transmission_zeros_diagonal_two_state() {
    // G(s) = 3/(s+1) - 1/(s+3) = 2(s+4) / ((s+1)(s+3))
    // Zero at s = -4, verifiable by hand:
    //   G(s) = C(sI-A)^{-1}B with A = diag(-1,-3), B = [1;1], C = [3,-1]
    //   numerator = 3(s+3) + (-1)(s+1) = 2s+8 = 2(s+4)
    caliburn::LinearSystem sys;
    sys.A = (Eigen::MatrixXd(2, 2) << -1, 0, 0, -3).finished();
    sys.B = (Eigen::MatrixXd(2, 1) <<  1, 1).finished();
    sys.C = (Eigen::MatrixXd(1, 2) <<  3, -1).finished();
    sys.D = Eigen::MatrixXd::Zero(1, 1);

    auto result = caliburn::computePoleZero(sys, 0, 0);

    ASSERT_EQ((int)result.zeros.size(), 1);
    ASSERT_NEAR(result.zeros[0].real(), -4.0, 1e-6);
    ASSERT_NEAR(result.zeros[0].imag(),  0.0, 1e-6);

    auto p1 = find_nearest(result.poles, {-1, 0});
    auto p2 = find_nearest(result.poles, {-3, 0});
    ASSERT_NEAR(p1.real(), -1.0, 1e-6);
    ASSERT_NEAR(p2.real(), -3.0, 1e-6);
}

int main() {
    test_first_order_pole();
    test_second_order_complex_poles();
    test_unstable_system();
    test_transmission_zeros();
    test_zeros_with_nonzero_D();
    test_root_locus_endpoints();
    test_state_feedback_locus();
    test_match_poles_uses_optimal_assignment();
    test_root_locus_sweep_branch_identity();
    test_transmission_zeros_diagonal_two_state();
    std::printf("All pole_zero tests passed.\n");
    return 0;
}
