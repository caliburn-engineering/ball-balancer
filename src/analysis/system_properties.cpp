// src/analysis/system_properties.cpp
#include "system_properties.h"
#include <Eigen/Eigenvalues>
#include <Eigen/QR>
#include <Eigen/SVD>
#include <complex>
#include <limits>

namespace caliburn {

namespace {

// Returns true if [ev*I - A, B] has full row rank (n), using JacobiSVD with a
// threshold of max(rows,cols) * epsilon * sigma_max -- the same rule MATLAB's
// rank() uses.
bool pbhFullRank(const LinearSystem& sys, std::complex<double> ev) {
    const int n = sys.states();
    const int m = sys.inputs();
    Eigen::MatrixXcd M(n, n + m);
    M.leftCols(n) = ev * Eigen::MatrixXcd::Identity(n, n)
                    - sys.A.cast<std::complex<double>>();
    M.rightCols(m) = sys.B.cast<std::complex<double>>();

    Eigen::JacobiSVD<Eigen::MatrixXcd> svd(
        M, Eigen::ComputeThinU | Eigen::ComputeThinV);
    if (svd.singularValues().size() == 0) return n == 0;
    const double tol = static_cast<double>(std::max(n, n + m))
                     * std::numeric_limits<double>::epsilon()
                     * svd.singularValues()(0);
    return static_cast<int>((svd.singularValues().array() > tol).count()) >= n;
}

}  // namespace

PropertyResult checkControllability(const LinearSystem& sys) {
    int n = sys.states();
    int m = sys.inputs();

    Eigen::MatrixXd ctrb(n, n * m);
    Eigen::MatrixXd Ak_B = sys.B;
    for (int k = 0; k < n; ++k) {
        ctrb.block(0, k * m, n, m) = Ak_B;
        Ak_B = sys.A * Ak_B;
    }

    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(ctrb);
    int rank = qr.rank();

    // PBH test for pass/fail: rank([λI - A, B]) == n for every eigenvalue λ.
    bool pass = false;
    Eigen::EigenSolver<Eigen::MatrixXd> es(sys.A, /*computeEigenvectors=*/false);
    if (es.info() == Eigen::Success) {
        pass = true;
        for (int i = 0; i < n && pass; ++i) {
            if (!pbhFullRank(sys, es.eigenvalues()(i))) pass = false;
        }
    }

    return {ctrb, rank, n, pass};
}

PropertyResult checkObservability(const LinearSystem& sys) {
    int n = sys.states();
    int p = sys.outputs();

    Eigen::MatrixXd obsv(n * p, n);
    Eigen::MatrixXd C_Ak = sys.C;
    for (int k = 0; k < n; ++k) {
        obsv.block(k * p, 0, p, n) = C_Ak;
        C_Ak = C_Ak * sys.A;
    }

    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(obsv);
    int rank = qr.rank();

    return {obsv, rank, n, rank == n};
}

bool isStabilizable(const LinearSystem& sys) {
    const int n = sys.states();
    if (n == 0) return true;

    Eigen::EigenSolver<Eigen::MatrixXd> es(sys.A, /*computeEigenvectors=*/false);
    if (es.info() != Eigen::Success) return false;

    // Only unstable (or marginally stable) modes need to be controllable.
    const double stab_tol = static_cast<double>(n)
                          * std::numeric_limits<double>::epsilon()
                          * std::max(1.0, sys.A.norm());

    for (int i = 0; i < n; ++i) {
        const std::complex<double> ev = es.eigenvalues()(i);
        if (ev.real() < -stab_tol) continue;  // stable mode -- no constraint
        if (!pbhFullRank(sys, ev)) return false;
    }
    return true;
}

}  // namespace caliburn
