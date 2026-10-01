// From caliburn reference/integrators/rk4.h @ c0b3246
#pragma once

#include <Eigen/Core>
#include <functional>

namespace caliburn {

/// Derivative function type: f(t, y) -> dy/dt
using DerivativeFn = std::function<Eigen::VectorXd(double, const Eigen::VectorXd&)>;

/// Advance state y at time t by one RK4 step of size h.
Eigen::VectorXd rk4_step(const Eigen::VectorXd& y, double t, double h,
                         const DerivativeFn& f);

/// Integrate from y0 at t0 for the given number of steps, returning the final state.
Eigen::VectorXd rk4_integrate(const Eigen::VectorXd& y0, double t0, double h,
                              int steps, const DerivativeFn& f);

}  // namespace caliburn
