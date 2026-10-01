// src/analysis/frequency_response.cpp
#include "frequency_response.h"
#include <Eigen/QR>
#include <cmath>

namespace caliburn {

std::complex<double> evalTransferFunction(
    const LinearSystem& sys, int output_i, int input_j,
    std::complex<double> s) {
    int n = sys.states();
    // n == 0 is a static gain.  ColPivHouseholderQR asserts on a 0x0 matrix in
    // a debug build; returning D directly is also exactly right.  See issue #5.
    if (n == 0) {
        return {sys.D(output_i, input_j), 0.0};
    }
    Eigen::MatrixXcd sI_A =
        s * Eigen::MatrixXcd::Identity(n, n) -
        sys.A.cast<std::complex<double>>();

    Eigen::VectorXcd b = sys.B.col(input_j).cast<std::complex<double>>();
    Eigen::VectorXcd x = sI_A.colPivHouseholderQr().solve(b);
    Eigen::VectorXcd c =
        sys.C.row(output_i).transpose().cast<std::complex<double>>();

    return c.dot(x) + std::complex<double>(sys.D(output_i, input_j), 0.0);
}

FrequencyResponse computeBode(
    const LinearSystem& sys, int output_i, int input_j,
    double freq_min_hz, double freq_max_hz, int num_points) {
    FrequencyResponse result;
    result.points.reserve(num_points);
    result.nyquist.reserve(num_points);

    double log_min = std::log10(freq_min_hz);
    double log_max = std::log10(freq_max_hz);

    double prev_phase_raw = 0.0;
    double phase_unwrap_offset = 0.0;

    for (int k = 0; k < num_points; ++k) {
        double log_f =
            log_min + (log_max - log_min) * k / (num_points - 1);
        double freq_hz = std::pow(10.0, log_f);
        double omega = 2.0 * M_PI * freq_hz;

        std::complex<double> s(0.0, omega);
        std::complex<double> G = evalTransferFunction(sys, output_i, input_j, s);

        double mag = std::abs(G);
        double magnitude_db = (mag > 0) ? 20.0 * std::log10(mag) : -300.0;
        double phase_raw = std::arg(G) * 180.0 / M_PI;

        // Phase unwrapping
        if (k > 0) {
            double diff = phase_raw - prev_phase_raw;
            if (diff > 180.0) phase_unwrap_offset -= 360.0;
            else if (diff < -180.0) phase_unwrap_offset += 360.0;
        }
        prev_phase_raw = phase_raw;
        double phase_deg = phase_raw + phase_unwrap_offset;

        result.points.push_back({freq_hz, magnitude_db, phase_deg});
        result.nyquist.push_back(G);
    }

    // Find gain crossover (|G| crosses 0 dB from above)
    for (int k = 1; k < num_points; ++k) {
        double m0 = result.points[k - 1].magnitude_db;
        double m1 = result.points[k].magnitude_db;
        if (m0 >= 0.0 && m1 < 0.0) {
            double t = -m0 / (m1 - m0);
            double lf0 = std::log10(result.points[k - 1].freq_hz);
            double lf1 = std::log10(result.points[k].freq_hz);
            result.gain_crossover_hz = std::pow(10.0, lf0 + t * (lf1 - lf0));
            double p0 = result.points[k - 1].phase_deg;
            double p1 = result.points[k].phase_deg;
            result.phase_margin_deg = (p0 + t * (p1 - p0)) + 180.0;
            break;
        }
    }

    // Find phase crossover (phase crosses -180° from above)
    for (int k = 1; k < num_points; ++k) {
        double p0 = result.points[k - 1].phase_deg;
        double p1 = result.points[k].phase_deg;
        if (p0 >= -180.0 && p1 < -180.0) {
            double t = (-180.0 - p0) / (p1 - p0);
            double lf0 = std::log10(result.points[k - 1].freq_hz);
            double lf1 = std::log10(result.points[k].freq_hz);
            result.phase_crossover_hz = std::pow(10.0, lf0 + t * (lf1 - lf0));
            double m0 = result.points[k - 1].magnitude_db;
            double m1 = result.points[k].magnitude_db;
            result.gain_margin_db = -(m0 + t * (m1 - m0));
            break;
        }
    }

    return result;
}

int nyquistEncirclements(const std::vector<std::complex<double>>& nyquist) {
    if (static_cast<int>(nyquist.size()) < 2) return 0;
    int n = static_cast<int>(nyquist.size());

    // Build the full closed Nyquist contour and compute the winding number of
    // G(jω) − (−1) around the origin.  The contour is:
    //   1. Positive frequencies:  nyquist[0] … nyquist[n-1]       (ω: 0 → +∞)
    //   2. Transition at ω → ∞:  nyquist[n-1] → conj(nyquist[n-1])
    //   3. Negative frequencies:  conj(nyquist[n-1]) … conj(nyquist[0])
    //   4. Closure at ω → 0:     conj(nyquist[0]) → nyquist[0]
    // For a real-coefficient, strictly proper system the DC gain is real so
    // step 4 contributes zero; step 2 is negligible when G(j∞) ≈ 0.

    const std::complex<double> critical(-1.0, 0.0);

    auto unwrapped_diff = [](double prev_arg, double curr_arg) -> double {
        double d = curr_arg - prev_arg;
        while (d > M_PI) d -= 2.0 * M_PI;
        while (d < -M_PI) d += 2.0 * M_PI;
        return d;
    };

    double total = 0.0;
    double prev = std::arg(nyquist[0] - critical);

    for (int k = 1; k < n; ++k) {
        double curr = std::arg(nyquist[k] - critical);
        total += unwrapped_diff(prev, curr);
        prev = curr;
    }

    {
        double curr = std::arg(std::conj(nyquist[n - 1]) - critical);
        total += unwrapped_diff(prev, curr);
        prev = curr;
    }

    for (int k = n - 2; k >= 0; --k) {
        double curr = std::arg(std::conj(nyquist[k]) - critical);
        total += unwrapped_diff(prev, curr);
        prev = curr;
    }

    {
        double curr = std::arg(nyquist[0] - critical);
        total += unwrapped_diff(prev, curr);
    }

    // CCW positive convention → flip sign for CW (Nyquist criterion's N).
    return -static_cast<int>(std::round(total / (2.0 * M_PI)));
}

}  // namespace caliburn
