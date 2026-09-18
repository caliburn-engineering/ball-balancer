// src/analysis/time_response.h
#pragma once

#include "linear_system.h"
#include <Eigen/Core>
#include <vector>

namespace caliburn {

struct TimePoint {
    double time;
    Eigen::VectorXd state;
    Eigen::VectorXd output;
    Eigen::VectorXd input;
};

struct TimeResponse {
    std::vector<TimePoint> points;
};

TimeResponse computeStepResponse(
    const LinearSystem& sys, int input_j,
    double amplitude, double duration, double dt);

TimeResponse computeImpulseResponse(
    const LinearSystem& sys, int input_j,
    double amplitude, double duration, double dt);

TimeResponse computeRampResponse(
    const LinearSystem& sys, int input_j,
    double slope, double duration, double dt);

}  // namespace caliburn
