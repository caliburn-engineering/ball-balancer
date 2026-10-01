// From caliburn reference/models/linearizer.h @ 2f7a6d0
#pragma once

#include "linear_system.h"

#include <Eigen/Core>
#include <functional>

namespace caliburn {

/// Nonlinear system: f(state, input) -> state_derivative
using NonlinearFn = std::function<Eigen::VectorXd(
    const Eigen::VectorXd& state,
    const Eigen::VectorXd& input)>;

/// Compute A = df/dx, B = df/du at (x0, u0) via central finite differences.
/// Returns LinearSystem with C = I, D = 0 (full state output).
LinearSystem linearize(
    const NonlinearFn& f,
    const Eigen::VectorXd& x0,
    const Eigen::VectorXd& u0,
    double epsilon = 1e-6);

/// Overload: caller provides C and D.
LinearSystem linearize(
    const NonlinearFn& f,
    const Eigen::VectorXd& x0,
    const Eigen::VectorXd& u0,
    const Eigen::MatrixXd& C,
    const Eigen::MatrixXd& D,
    double epsilon = 1e-6);

/// Element-wise comparison of analytical vs numerical linearization.
struct ValidationResult {
    Eigen::MatrixXd A_error;  // |analytical.A - numerical.A|
    Eigen::MatrixXd B_error;  // |analytical.B - numerical.B|
    double max_A_error;
    double max_B_error;
    bool pass;                // both maxes < tolerance
};

ValidationResult validate(
    const LinearSystem& analytical,
    const NonlinearFn& f,
    const Eigen::VectorXd& x0,
    const Eigen::VectorXd& u0,
    double tol = 1e-4);

}  // namespace caliburn
