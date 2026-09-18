#include "comparison_sim.h"
#include "rk4.h"

#include <cmath>

namespace caliburn {

Eigen::VectorXd evaluateInput(const InputSignal& signal, double t, int num_inputs) {
    return std::visit([&](auto&& sig) -> Eigen::VectorXd {
        using T = std::decay_t<decltype(sig)>;
        if constexpr (std::is_same_v<T, StepInput>) {
            if (t >= sig.start_time)
                return sig.amplitude;
            return Eigen::VectorXd::Zero(num_inputs);
        } else if constexpr (std::is_same_v<T, RampInput>) {
            if (t >= sig.start_time)
                return sig.slope * (t - sig.start_time);
            return Eigen::VectorXd::Zero(num_inputs);
        } else {  // CustomInput
            return sig(t);
        }
    }, signal);
}

ComparisonResult compareModels(
    const NonlinearFn& f_nl,
    const LinearSystem& lin,
    const Eigen::VectorXd& x0,
    const InputSignal& input,
    double dt,
    double duration)
{
    int num_steps = static_cast<int>(std::ceil(duration / dt));
    int n = static_cast<int>(x0.size());
    int m = lin.inputs();

    ComparisonResult result;
    result.time.reserve(num_steps + 1);
    result.nonlinear_states.reserve(num_steps + 1);
    result.linear_states.reserve(num_steps + 1);
    result.inputs.reserve(num_steps + 1);

    Eigen::VectorXd x_nl = x0;
    Eigen::VectorXd x_lin = x0;

    for (int i = 0; i <= num_steps; ++i) {
        double t = i * dt;
        Eigen::VectorXd u = evaluateInput(input, t, m);

        result.time.push_back(t);
        result.nonlinear_states.push_back(x_nl);
        result.linear_states.push_back(x_lin);
        result.inputs.push_back(u);

        if (i < num_steps) {
            // RK4 step for nonlinear model
            DerivativeFn f_nl_rk4 = [&](double /*t*/, const Eigen::VectorXd& x) {
                return f_nl(x, u);
            };
            x_nl = rk4_step(x_nl, t, dt, f_nl_rk4);

            // RK4 step for linear model: x_dot = A*x + B*u
            DerivativeFn f_lin_rk4 = [&](double /*t*/, const Eigen::VectorXd& x) {
                return Eigen::VectorXd(lin.A * x + lin.B * u);
            };
            x_lin = rk4_step(x_lin, t, dt, f_lin_rk4);
        }
    }

    return result;
}

}  // namespace caliburn
