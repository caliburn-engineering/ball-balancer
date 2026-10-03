// tests/test_frequency_response.cpp
#include "test_helpers.h"
#include "analysis/frequency_response.h"
#include <cmath>
#include <cstdio>

// First-order system: G(s) = 1/(s+1)
// A=[-1], B=[1], C=[1], D=[0]
caliburn::LinearSystem make_first_order() {
    caliburn::LinearSystem sys;
    sys.A = (Eigen::MatrixXd(1, 1) << -1).finished();
    sys.B = (Eigen::MatrixXd(1, 1) << 1).finished();
    sys.C = (Eigen::MatrixXd(1, 1) << 1).finished();
    sys.D = Eigen::MatrixXd::Zero(1, 1);
    return sys;
}

// G(s) = (s+3)/((s+1)(s+2)) — has a zero at s=-3
// A=[0 1; -2 -3], B=[0; 1], C=[3 1], D=[0]
caliburn::LinearSystem make_system_with_zero() {
    caliburn::LinearSystem sys;
    sys.A = (Eigen::MatrixXd(2, 2) << 0, 1, -2, -3).finished();
    sys.B = (Eigen::MatrixXd(2, 1) << 0, 1).finished();
    sys.C = (Eigen::MatrixXd(1, 2) << 3, 1).finished();
    sys.D = Eigen::MatrixXd::Zero(1, 1);
    return sys;
}

void test_eval_dc_gain() {
    auto sys = make_first_order();
    std::complex<double> s(0.0, 0.0);
    auto G = caliburn::evalTransferFunction(sys, 0, 0, s);
    ASSERT_NEAR(G.real(), 1.0, 1e-10);
    ASSERT_NEAR(G.imag(), 0.0, 1e-10);
}

void test_eval_at_bandwidth() {
    auto sys = make_first_order();
    std::complex<double> s(0.0, 1.0);
    auto G = caliburn::evalTransferFunction(sys, 0, 0, s);
    ASSERT_NEAR(std::abs(G), 1.0 / std::sqrt(2.0), 1e-10);
    ASSERT_NEAR(std::arg(G) * 180.0 / M_PI, -45.0, 1e-6);
}

void test_bode_first_order() {
    auto sys = make_first_order();
    auto resp = caliburn::computeBode(sys, 0, 0, 0.01, 100.0, 500);
    ASSERT_TRUE(resp.points.size() == 500);
    ASSERT_NEAR(resp.points[0].magnitude_db, 0.0, 0.1);

    double target_f = 1.0 / (2.0 * M_PI);
    int best_idx = 0;
    double best_diff = 1e10;
    for (int k = 0; k < (int)resp.points.size(); ++k) {
        double diff = std::abs(resp.points[k].freq_hz - target_f);
        if (diff < best_diff) {
            best_diff = diff;
            best_idx = k;
        }
    }
    ASSERT_NEAR(resp.points[best_idx].magnitude_db, -3.01, 0.5);
    ASSERT_NEAR(resp.points[best_idx].phase_deg, -45.0, 2.0);

    ASSERT_TRUE(resp.points.back().magnitude_db < -30.0);
    ASSERT_NEAR(resp.points.back().phase_deg, -90.0, 5.0);
}

void test_bode_nyquist_data() {
    auto sys = make_first_order();
    auto resp = caliburn::computeBode(sys, 0, 0, 0.01, 100.0, 100);
    ASSERT_EQ((int)resp.nyquist.size(), 100);
    ASSERT_NEAR(resp.nyquist[0].real(), 1.0, 0.1);
    ASSERT_NEAR(resp.nyquist[0].imag(), 0.0, 0.1);
}

void test_phase_continuity() {
    auto sys = make_system_with_zero();
    auto resp = caliburn::computeBode(sys, 0, 0, 0.01, 100.0, 500);
    for (int k = 1; k < (int)resp.points.size(); ++k) {
        double diff = std::abs(resp.points[k].phase_deg - resp.points[k - 1].phase_deg);
        ASSERT_TRUE(diff < 30.0);
    }
}

void test_margins_stable_system() {
    // G(s) = 100/((s+1)(s+5)(s+10)) — stable, gain > 0 dB at low freq, finite margins
    // DC gain = 100/50 = 2 → ~6 dB (crosses 0 dB)
    // CCF: den = s³ + 16s² + 65s + 50
    // A = [0 1 0; 0 0 1; -50 -65 -16], B = [0; 0; 1], C = [100 0 0], D = [0]
    caliburn::LinearSystem sys;
    sys.A = (Eigen::MatrixXd(3, 3) << 0, 1, 0, 0, 0, 1, -50, -65, -16)
                .finished();
    sys.B = (Eigen::MatrixXd(3, 1) << 0, 0, 1).finished();
    sys.C = (Eigen::MatrixXd(1, 3) << 100, 0, 0).finished();
    sys.D = Eigen::MatrixXd::Zero(1, 1);

    auto resp = caliburn::computeBode(sys, 0, 0, 0.001, 100.0, 5000);

    ASSERT_NEAR(resp.points[0].magnitude_db, 6.02, 0.5);
    ASSERT_TRUE(resp.gain_crossover_hz > 0);
    ASSERT_TRUE(resp.phase_margin_deg > 0);
    ASSERT_TRUE(resp.phase_margin_deg < 180);
    ASSERT_TRUE(resp.phase_crossover_hz > 0);
    ASSERT_TRUE(resp.gain_margin_db > 0);
}

void test_margins_type1_system() {
    // G(s) = 1/(s(s+1)(s+2)) — type-1 system
    // CCF: den = s³ + 3s² + 2s
    // A = [0 1 0; 0 0 1; 0 -2 -3], B = [0; 0; 1], C = [1 0 0], D = [0]
    caliburn::LinearSystem sys;
    sys.A = (Eigen::MatrixXd(3, 3) << 0, 1, 0, 0, 0, 1, 0, -2, -3).finished();
    sys.B = (Eigen::MatrixXd(3, 1) << 0, 0, 1).finished();
    sys.C = (Eigen::MatrixXd(1, 3) << 1, 0, 0).finished();
    sys.D = Eigen::MatrixXd::Zero(1, 1);

    auto resp = caliburn::computeBode(sys, 0, 0, 0.001, 100.0, 5000);

    ASSERT_NEAR(resp.gain_margin_db, 15.56, 1.0);
    ASSERT_NEAR(resp.phase_crossover_hz, std::sqrt(2.0) / (2.0 * M_PI), 0.02);
    ASSERT_TRUE(resp.gain_crossover_hz > 0);
    ASSERT_TRUE(resp.phase_margin_deg > 0);
    ASSERT_TRUE(resp.phase_margin_deg < 90);
}

// G(s) = 1/((s+1)(s+2)) in CCF: denomintor s²+3s+2.
// A=[0 1;-2 -3], B=[0;1], C=[1 0], D=[0].
caliburn::LinearSystem make_second_order() {
    caliburn::LinearSystem sys;
    sys.A = (Eigen::MatrixXd(2, 2) << 0, 1, -2, -3).finished();
    sys.B = (Eigen::MatrixXd(2, 1) << 0, 1).finished();
    sys.C = (Eigen::MatrixXd(1, 2) << 1, 0).finished();
    sys.D = Eigen::MatrixXd::Zero(1, 1);
    return sys;
}

// Verify computeBode magnitude at ω=1 rad/s against the analytic value.
// G(s) = 1/((s+1)(s+2)) → |G(j1)| = 1/√((1+1)(4+1)) = 1/√10.
void test_bode_analytic_second_order() {
    auto sys = make_second_order();
    // Verify via evalTransferFunction at jω = j*1 rad/s.
    std::complex<double> s(0.0, 1.0);
    auto G = caliburn::evalTransferFunction(sys, 0, 0, s);
    ASSERT_NEAR(std::abs(G), 1.0 / std::sqrt(10.0), 1e-10);
    ASSERT_NEAR(20.0 * std::log10(std::abs(G)), -10.0, 1e-6);

    // Verify the same point is reproduced by computeBode (nearest grid point).
    auto resp = caliburn::computeBode(sys, 0, 0, 0.001, 100.0, 2000);
    double target_f = 1.0 / (2.0 * M_PI);
    int best = 0;
    double best_diff = 1e10;
    for (int k = 0; k < (int)resp.points.size(); ++k) {
        double d = std::abs(resp.points[k].freq_hz - target_f);
        if (d < best_diff) { best_diff = d; best = k; }
    }
    ASSERT_NEAR(resp.points[best].magnitude_db, -10.0, 0.5);
}

// G(s) = 1/((s+1)(s+2)): stable OL (P=0), stable CL (Z=0) → N=0.
void test_nyquist_stable_no_encirclement() {
    auto sys = make_second_order();
    auto resp = caliburn::computeBode(sys, 0, 0, 0.001, 10000.0, 5000);
    ASSERT_EQ(caliburn::nyquistEncirclements(resp.nyquist), 0);
}

// G(s) = 100/((s+1)(s+2)(s+3)): P=0, gain margin negative.
// CL characteristic equation s³+6s²+11s+106 — Routh gives Z=2.
// Nyquist criterion: N = Z−P = 2.
void test_nyquist_two_cw_encirclements() {
    // CCF for 100/(s³+6s²+11s+6): A=[0 1 0;0 0 1;-6 -11 -6], B=[0;0;1], C=[100 0 0].
    caliburn::LinearSystem sys;
    sys.A = (Eigen::MatrixXd(3, 3) << 0, 1, 0, 0, 0, 1, -6, -11, -6).finished();
    sys.B = (Eigen::MatrixXd(3, 1) << 0, 0, 1).finished();
    sys.C = (Eigen::MatrixXd(1, 3) << 100, 0, 0).finished();
    sys.D = Eigen::MatrixXd::Zero(1, 1);
    auto resp = caliburn::computeBode(sys, 0, 0, 0.001, 10000.0, 5000);
    ASSERT_EQ(caliburn::nyquistEncirclements(resp.nyquist), 2);
}

// G(s) = 10/(s−1): one unstable OL pole (P=1).
// CL pole at s=−9 (stable), so Z=0 → N = Z−P = −1 (one CCW encirclement).
void test_nyquist_one_ccw_encirclement_unstable_ol() {
    // A=[1], B=[1], C=[10], D=[0].
    caliburn::LinearSystem sys;
    sys.A = (Eigen::MatrixXd(1, 1) << 1).finished();
    sys.B = (Eigen::MatrixXd(1, 1) << 1).finished();
    sys.C = (Eigen::MatrixXd(1, 1) << 10).finished();
    sys.D = Eigen::MatrixXd::Zero(1, 1);
    auto resp = caliburn::computeBode(sys, 0, 0, 0.001, 10000.0, 5000);
    ASSERT_EQ(caliburn::nyquistEncirclements(resp.nyquist), -1);
}

int main() {
    test_eval_dc_gain();
    test_eval_at_bandwidth();
    test_bode_first_order();
    test_bode_nyquist_data();
    test_phase_continuity();
    test_margins_stable_system();
    test_margins_type1_system();
    test_bode_analytic_second_order();
    test_nyquist_stable_no_encirclement();
    test_nyquist_two_cw_encirclements();
    test_nyquist_one_ccw_encirclement_unstable_ol();
    std::printf("All frequency_response tests passed.\n");
    return 0;
}
