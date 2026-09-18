// src/analysis/loop_diagnostics.h
#pragma once

#include "linear_system.h"
#include "controller_builders.h"

#include <Eigen/Core>
#include <string>
#include <vector>

namespace caliburn {

// The two readouts are mutually exclusive: a square paired sub-matrix of size
// >= 2 gives the RGA, everything else gives Channel Share.  One panel slot.
struct LoopDiagnostics {
    // --- RGA (pairing diagnostic) ---
    bool has_rga = false;
    Eigen::MatrixXd lambda;            // q x q, the real part of the RGA
    std::vector<int> sub_out, sub_in;  // sub-matrix index maps into (p, m)
    double rga_number = 0.0;           // ||L - Pi||_sum against the actual pairing

    // --- Channel Share (NOT a pairing measure) ---
    bool has_share = false;
    int share_in = 0;                  // the paired input the shares are of
    std::vector<double> share;         // per plant output, sums to 1
    std::vector<double> mag_db;        // |g_i(jw)| in dB, UNSCALED, per output

    // --- per loop ---
    std::vector<double> loop_lambda;   // NaN where there is no RGA
    std::vector<char> loop_dead;       // structurally dead paired channel

    bool scales_at_default = true;
    std::string headline;
};

// G(jw) as a complex matrix, one evalTransferFunction per channel.
Eigen::MatrixXcd evalPlantMatrix(const LinearSystem& sys, double freq_hz);

// |g_ij(jw)| below tol at EVERY w on the Bode grid.  Positive rescaling cannot
// turn a zero into a nonzero at any frequency, so this test is immune to the
// scale-dependence of Channel Share — unlike any threshold on the share itself,
// or a ratio against the column maximum.
bool isStructurallyDeadChannel(const LinearSystem& sys, int output_i, int input_j,
                               double freq_min_hz, double freq_max_hz,
                               int num_points, double tol = 1e-9);

LoopDiagnostics computeLoopDiagnostics(const LinearSystem& plant,
                                       const std::vector<Loop>& loops,
                                       const std::vector<float>& output_scales,
                                       double freq_hz,
                                       double freq_min_hz, double freq_max_hz,
                                       int num_points);

// The RGA number alone at one w, for the sweep strip.  NaN where no RGA exists.
double rgaNumberAt(const LinearSystem& plant, const std::vector<Loop>& loops,
                   double freq_hz);

}  // namespace caliburn
