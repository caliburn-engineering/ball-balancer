// src/analysis/system_properties.h
#pragma once

#include "linear_system.h"
#include <Eigen/Core>

namespace caliburn {

struct PropertyResult {
    Eigen::MatrixXd matrix;  // the controllability or observability matrix
    int rank;
    int required_rank;       // = n (number of states)
    bool pass;               // rank == required_rank
};

// Controllability matrix: [B, AB, A²B, ..., A^(n-1)B]
PropertyResult checkControllability(const LinearSystem& sys);

// Observability matrix: [C; CA; CA²; ...; CA^(n-1)]
PropertyResult checkObservability(const LinearSystem& sys);

}  // namespace caliburn
