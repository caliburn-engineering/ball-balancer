// src/analysis/time_response.cpp
#include "time_response.h"
#include "rk4.h"

namespace caliburn {

TimeResponse computeStepResponse(
    const LinearSystem& sys, int input_j,
    double amplitude, double duration, double dt) {
    int n = sys.states();
    int m = sys.inputs();
    TimeResponse result;

    Eigen::VectorXd x = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd u = Eigen::VectorXd::Zero(m);
    u(input_j) = amplitude;

    int steps = static_cast<int>(duration / dt);
    result.points.reserve(steps + 1);

    for (int k = 0; k <= steps; ++k) {
        double t = k * dt;
        Eigen::VectorXd y = sys.C * x + sys.D * u;
        result.points.push_back({t, x, y, u});

        if (k < steps) {
            DerivativeFn deriv = [&](double, const Eigen::VectorXd& state) -> Eigen::VectorXd {
                return sys.A * state + sys.B * u;
            };
            x = rk4_step(x, t, dt, deriv);
        }
    }
    return result;
}

TimeResponse computeImpulseResponse(
    const LinearSystem& sys, int input_j,
    double amplitude, double duration, double dt) {
    int n = sys.states();
    int m = sys.inputs();
    TimeResponse result;

    Eigen::VectorXd x = Eigen::VectorXd::Zero(n);
    int steps = static_cast<int>(duration / dt);
    result.points.reserve(steps + 1);

    for (int k = 0; k <= steps; ++k) {
        double t = k * dt;
        Eigen::VectorXd u = Eigen::VectorXd::Zero(m);
        if (k == 0) u(input_j) = amplitude / dt;

        Eigen::VectorXd y = sys.C * x + sys.D * u;
        result.points.push_back({t, x, y, u});

        if (k < steps) {
            DerivativeFn deriv = [&](double, const Eigen::VectorXd& state) -> Eigen::VectorXd {
                return sys.A * state + sys.B * u;
            };
            x = rk4_step(x, t, dt, deriv);
        }
    }
    return result;
}

TimeResponse computeRampResponse(
    const LinearSystem& sys, int input_j,
    double slope, double duration, double dt) {
    int n = sys.states();
    int m = sys.inputs();
    TimeResponse result;

    Eigen::VectorXd x = Eigen::VectorXd::Zero(n);
    int steps = static_cast<int>(duration / dt);
    result.points.reserve(steps + 1);

    for (int k = 0; k <= steps; ++k) {
        double t = k * dt;
        Eigen::VectorXd u = Eigen::VectorXd::Zero(m);
        u(input_j) = slope * t;

        Eigen::VectorXd y = sys.C * x + sys.D * u;
        result.points.push_back({t, x, y, u});

        if (k < steps) {
            DerivativeFn deriv = [&](double t_inner,
                                     const Eigen::VectorXd& state) -> Eigen::VectorXd {
                Eigen::VectorXd u_inner = Eigen::VectorXd::Zero(m);
                u_inner(input_j) = slope * t_inner;
                return sys.A * state + sys.B * u_inner;
            };
            x = rk4_step(x, t, dt, deriv);
        }
    }
    return result;
}

}  // namespace caliburn
