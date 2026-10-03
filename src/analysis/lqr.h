// src/analysis/lqr.h
#pragma once

#include "linear_system.h"
#include <Eigen/Core>
#include <complex>
#include <string>
#include <vector>

namespace caliburn {

struct LqrResult {
    Eigen::MatrixXd K;  // m x n optimal state-feedback gain, u = -Kx
    Eigen::MatrixXd P;  // n x n stabilizing solution of the CARE
    std::vector<std::complex<double>> closed_loop_poles;  // eig(A - BK)
    // Relative Riccati residual after Newton/Kleinman refinement:
    //   ||A'P + PA - PBR^-1B'P + Q|| / (||A'P|| + ||PA|| + ||PBR^-1B'P|| + ||Q||)
    double residual = 0.0;
    // Residual of the sign-function solution, before refinement.  On the
    // record so the test can report it without needing to re-run the solver.
    double pre_refinement_residual = 0.0;
    bool success = false;
    std::string error;
};

// Continuous-time LQR: minimize J = integral(x'Qx + u'Ru) dt subject to
// x' = Ax + Bu.  Returns K = R^-1 B' P, where P solves the continuous
// algebraic Riccati equation
//
//     A'P + PA - PBR^-1B'P + Q = 0
//
// via the matrix sign function.  Only `sys.A` and `sys.B` are read; C and D
// play no part in a regulator.
//
// On failure `success` is false, `error` states why, and K, P and the pole
// list are empty.  Failure is a rejected input, a non-converging iteration, or
// a relative Riccati residual above `residual_gate` after refinement — never a
// silently wrong gain.
//
// Normal solutions land at 1e-12 or below, so the default gate trips only on
// genuine breakage.  Callers leave it alone; a test lowers it to make the gate
// fire on a well-conditioned plant.
constexpr double kLqrResidualGate = 1e-8;

LqrResult computeLQR(const LinearSystem& sys,
                     const Eigen::MatrixXd& Q,
                     const Eigen::MatrixXd& R,
                     double residual_gate = kLqrResidualGate);

}  // namespace caliburn
