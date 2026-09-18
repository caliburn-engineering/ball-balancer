// src/analysis/system_connect.h
#pragma once

#include "linear_system.h"
#include <Eigen/Core>

namespace caliburn {

LinearSystem seriesConnect(const LinearSystem& sys1, const LinearSystem& sys2);
LinearSystem feedbackConnect(const LinearSystem& open_loop);
LinearSystem stateFeedbackClose(
    const LinearSystem& plant, const Eigen::MatrixXd& K);

}  // namespace caliburn
