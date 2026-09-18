#pragma once

#include "linear_system.h"
#include "linearizer.h"

#include <Eigen/Core>
#include <functional>
#include <variant>
#include <vector>

namespace caliburn {

// --- Input signal types ---

struct StepInput {
    Eigen::VectorXd amplitude;  // step size per input channel
    double start_time = 0.0;
};

struct RampInput {
    Eigen::VectorXd slope;  // slope per input channel
    double start_time = 0.0;
};

using CustomInput = std::function<Eigen::VectorXd(double t)>;
using InputSignal = std::variant<StepInput, RampInput, CustomInput>;

/// Evaluate an InputSignal at time t
Eigen::VectorXd evaluateInput(const InputSignal& signal, double t, int num_inputs);

// --- Comparison result ---

struct ComparisonResult {
    std::vector<double> time;
    std::vector<Eigen::VectorXd> nonlinear_states;
    std::vector<Eigen::VectorXd> linear_states;
    std::vector<Eigen::VectorXd> inputs;
};

/// Run both nonlinear and linear models with the same input and initial condition.
/// Both use RK4 integration with the same dt.
ComparisonResult compareModels(
    const NonlinearFn& f_nl,
    const LinearSystem& lin,
    const Eigen::VectorXd& x0,
    const InputSignal& input,
    double dt,
    double duration);

}  // namespace caliburn
