// src/analysis/frequency_response.h
#pragma once

#include "linear_system.h"
#include <Eigen/Core>
#include <complex>
#include <limits>
#include <vector>

namespace caliburn {

struct BodePoint {
    double freq_hz;
    double magnitude_db;
    double phase_deg;
};

struct FrequencyResponse {
    std::vector<BodePoint> points;
    std::vector<std::complex<double>> nyquist;  // raw G(jω) at each freq point
    double gain_margin_db = std::numeric_limits<double>::infinity();
    double phase_margin_deg = std::numeric_limits<double>::infinity();
    double gain_crossover_hz = -1.0;   // where |G| = 0 dB
    double phase_crossover_hz = -1.0;  // where ∠G = -180°
};

// Evaluate single SISO channel (input j → output i) at complex frequency s.
// Computes G_ij(s) = C_i (sI - A)⁻¹ B_j + D_ij.
std::complex<double> evalTransferFunction(
    const LinearSystem& sys, int output_i, int input_j,
    std::complex<double> s);

// Evaluate channel (i, j) over a log-spaced frequency grid.
// Phase is continuously unwrapped. Gain/phase margins computed.
FrequencyResponse computeBode(
    const LinearSystem& sys, int output_i, int input_j,
    double freq_min_hz, double freq_max_hz, int num_points);

}  // namespace caliburn
