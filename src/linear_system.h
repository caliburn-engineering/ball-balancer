// From caliburn reference/models/linear_system.h @ b457f57
#pragma once

#include <Eigen/Core>

namespace caliburn {

struct LinearSystem {
    Eigen::MatrixXd A;  // n x n state matrix
    Eigen::MatrixXd B;  // n x m input matrix
    Eigen::MatrixXd C;  // p x n output matrix
    Eigen::MatrixXd D;  // p x m feedthrough matrix

    int states() const { return static_cast<int>(A.rows()); }
    int inputs() const { return static_cast<int>(B.cols()); }
    int outputs() const { return static_cast<int>(C.rows()); }
};

}  // namespace caliburn
