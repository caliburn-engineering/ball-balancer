// From caliburn reference/integrators/rk4.cpp @ c0b3246
#include "rk4.h"

namespace caliburn {

Eigen::VectorXd rk4_step(const Eigen::VectorXd& y, double t, double h,
                         const DerivativeFn& f) {
    const Eigen::VectorXd k1 = f(t, y);
    const Eigen::VectorXd k2 = f(t + h / 2.0, y + (h / 2.0) * k1);
    const Eigen::VectorXd k3 = f(t + h / 2.0, y + (h / 2.0) * k2);
    const Eigen::VectorXd k4 = f(t + h, y + h * k3);

    return y + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
}

Eigen::VectorXd rk4_integrate(const Eigen::VectorXd& y0, double t0, double h,
                              int steps, const DerivativeFn& f) {
    Eigen::VectorXd y = y0;
    double t = t0;

    for (int i = 0; i < steps; ++i) {
        y = rk4_step(y, t, h, f);
        t += h;
    }

    return y;
}

}  // namespace caliburn
