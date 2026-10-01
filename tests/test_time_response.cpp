// tests/test_time_response.cpp
#include "test_helpers.h"
#include "analysis/time_response.h"
#include "assert_rel.h"
#include <cmath>
#include <cstdio>

caliburn::LinearSystem make_first_order_sys() {
    caliburn::LinearSystem sys;
    sys.A = (Eigen::MatrixXd(1, 1) << -1).finished();
    sys.B = (Eigen::MatrixXd(1, 1) << 1).finished();
    sys.C = (Eigen::MatrixXd(1, 1) << 1).finished();
    sys.D = Eigen::MatrixXd::Zero(1, 1);
    return sys;
}

caliburn::LinearSystem make_integrator() {
    caliburn::LinearSystem sys;
    sys.A = Eigen::MatrixXd::Zero(1, 1);
    sys.B = (Eigen::MatrixXd(1, 1) << 1).finished();
    sys.C = (Eigen::MatrixXd(1, 1) << 1).finished();
    sys.D = Eigen::MatrixXd::Zero(1, 1);
    return sys;
}

void test_step_first_order() {
    auto sys = make_first_order_sys();
    auto resp = caliburn::computeStepResponse(sys, 0, 1.0, 5.0, 0.001);
    ASSERT_TRUE(resp.points.size() > 100);
    ASSERT_NEAR(resp.points[0].output(0), 0.0, 1e-6);
    int idx_1s = static_cast<int>(1.0 / 0.001);
    ASSERT_NEAR(resp.points[idx_1s].output(0), 1.0 - std::exp(-1.0), 1e-3);
    ASSERT_NEAR(resp.points.back().output(0), 1.0, 0.01);
}

void test_step_amplitude() {
    auto sys = make_first_order_sys();
    auto resp = caliburn::computeStepResponse(sys, 0, 3.0, 10.0, 0.001);
    ASSERT_NEAR(resp.points.back().output(0), 3.0, 0.01);
}

void test_impulse_first_order() {
    auto sys = make_first_order_sys();
    auto resp = caliburn::computeImpulseResponse(sys, 0, 1.0, 5.0, 0.001);
    int idx_1s = static_cast<int>(1.0 / 0.001);
    ASSERT_NEAR(resp.points[idx_1s].output(0), std::exp(-1.0), 0.01);
    ASSERT_NEAR(resp.points.back().output(0), 0.0, 0.01);
}

void test_step_integrator() {
    auto sys = make_integrator();
    auto resp = caliburn::computeStepResponse(sys, 0, 1.0, 5.0, 0.001);
    int idx_2s = static_cast<int>(2.0 / 0.001);
    ASSERT_NEAR(resp.points[idx_2s].output(0), 2.0, 0.01);
    ASSERT_NEAR(resp.points.back().output(0), 5.0, 0.01);
}

void test_ramp_integrator() {
    auto sys = make_integrator();
    auto resp = caliburn::computeRampResponse(sys, 0, 1.0, 5.0, 0.001);
    int idx_2s = static_cast<int>(2.0 / 0.001);
    ASSERT_NEAR(resp.points[idx_2s].output(0), 2.0, 0.02);
    int idx_4s = static_cast<int>(4.0 / 0.001);
    ASSERT_NEAR(resp.points[idx_4s].output(0), 8.0, 0.05);
}

void test_time_response_dimensions() {
    caliburn::LinearSystem sys;
    sys.A = Eigen::MatrixXd::Zero(2, 2);
    sys.A(0, 0) = -1;
    sys.A(1, 1) = -2;
    sys.B = Eigen::MatrixXd::Identity(2, 2);
    sys.C = Eigen::MatrixXd::Identity(2, 2);
    sys.D = Eigen::MatrixXd::Zero(2, 2);

    auto resp = caliburn::computeStepResponse(sys, 0, 1.0, 2.0, 0.01);
    ASSERT_EQ(resp.points[0].output.size(), 2);
    ASSERT_EQ(resp.points[0].input.size(), 2);
    ASSERT_EQ(resp.points[0].state.size(), 2);
    ASSERT_NEAR(resp.points[1].input(0), 1.0, 1e-10);
    ASSERT_NEAR(resp.points[1].input(1), 0.0, 1e-10);
}

// Decoupled 2-state system: A = diag(-1, -2), B = I, C = I, D = 0.
// Matrix exponential is expm(A*t) = diag(exp(-t), exp(-2t)).
caliburn::LinearSystem make_decoupled_second_order() {
    caliburn::LinearSystem sys;
    sys.A = (Eigen::MatrixXd(2, 2) << -1, 0, 0, -2).finished();
    sys.B = Eigen::MatrixXd::Identity(2, 2);
    sys.C = Eigen::MatrixXd::Identity(2, 2);
    sys.D = Eigen::MatrixXd::Zero(2, 2);
    return sys;
}

// Impulse j=0, amplitude a: x(0) = B[:,0]*a = (a,0).
// Exact: y(t) = expm(A*t) * x(0) = (a*exp(-t), 0).
void test_impulse_matches_matrix_exp() {
    auto sys = make_decoupled_second_order();
    double amplitude = 2.0;
    auto resp = caliburn::computeImpulseResponse(sys, 0, amplitude, 4.0, 0.001);
    for (int i = 0; i < (int)resp.points.size(); i += 500) {
        double t = resp.points[i].time;
        Eigen::VectorXd expected(2);
        expected(0) = amplitude * std::exp(-t);
        expected(1) = 0.0;
        ASSERT_MATRIX_REL_NEAR(resp.points[i].output, expected, 1e-3);
    }
}

// Step j=0, amplitude a: x(0)=0, u=(a,0).
// Exact: y_1(t) = a*(1-exp(-t)), y_2(t) = 0.
void test_step_matches_matrix_exp() {
    auto sys = make_decoupled_second_order();
    double amplitude = 3.0;
    auto resp = caliburn::computeStepResponse(sys, 0, amplitude, 4.0, 0.001);
    for (int i = 100; i < (int)resp.points.size(); i += 500) {
        double t = resp.points[i].time;
        Eigen::VectorXd expected(2);
        expected(0) = amplitude * (1.0 - std::exp(-t));
        expected(1) = 0.0;
        ASSERT_MATRIX_REL_NEAR(resp.points[i].output, expected, 1e-3);
    }
}

int main() {
    test_step_first_order();
    test_step_amplitude();
    test_impulse_first_order();
    test_step_integrator();
    test_ramp_integrator();
    test_time_response_dimensions();
    test_impulse_matches_matrix_exp();
    test_step_matches_matrix_exp();
    std::printf("All time_response tests passed.\n");
    return 0;
}
