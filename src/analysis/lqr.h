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
// list are empty.  Failure is a rejected input or a non-converging iteration,
// never a silently wrong gain.
LqrResult computeLQR(const LinearSystem& sys,
                     const Eigen::MatrixXd& Q,
                     const Eigen::MatrixXd& R);

}  // namespace caliburn
