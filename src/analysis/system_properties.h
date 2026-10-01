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

// Controllability matrix: [B, AB, A²B, ..., A^(n-1)B].
// `pass` uses the PBH/Hautus test (rank([λI-A, B]) == n for every eigenvalue
// λ of A), which is more reliable than Krylov rank at mixed units.  The matrix
// and rank fields are still the Krylov matrix and its rank, for display.
PropertyResult checkControllability(const LinearSystem& sys);

// Observability matrix: [C; CA; CA²; ...; CA^(n-1)]
PropertyResult checkObservability(const LinearSystem& sys);

// Stabilizability via PBH: returns true if every eigenvalue of A with
// Re(λ) >= 0 satisfies rank([λI - A, B]) == n.  A stabilizable-but-not-
// controllable plant has a stable uncontrollable mode and passes this test.
bool isStabilizable(const LinearSystem& sys);

}  // namespace caliburn
