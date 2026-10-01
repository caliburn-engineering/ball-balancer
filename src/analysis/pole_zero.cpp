// src/analysis/pole_zero.cpp
#include "pole_zero.h"
#include <Eigen/Eigenvalues>
#include <cmath>
#include <limits>
#include <vector>

namespace caliburn {

namespace {
// Hungarian algorithm: minimum total distance assignment of n current poles
// to n previous poles.  Returns perm where perm[i] = j means prev[i] is
// matched to poles[j].  O(n³), correct for any non-negative real costs.
//
// The greedy (nearest-neighbour) alternative is O(n²) but picks the locally
// nearest unmatched pole for each branch in order, which can assign a pole to
// the wrong branch when two branches pass near each other and the step is
// large: greedy "steals" the closer pole for an earlier branch, leaving a
// farther one for a later branch whose total cost would have been lower.
std::vector<int> hungarianAssign(
    const std::vector<std::complex<double>>& poles,
    const std::vector<std::complex<double>>& prev) {
    int n = static_cast<int>(poles.size());
    const double kInf = std::numeric_limits<double>::max() / 2.0;

    // u[i], v[j]: dual potentials (1-indexed; index 0 is a sentinel row).
    std::vector<double> u(n + 1, 0.0), v(n + 1, 0.0);
    // p[j]: which prev row is assigned to current column j (1-indexed).
    // way[j]: which column j0 caused column j's current tentative assignment.
    std::vector<int> p(n + 1, 0), way(n + 1, 0);

    for (int i = 1; i <= n; ++i) {
        p[0] = i;
        int j0 = 0;
        std::vector<double> minv(n + 1, kInf);
        std::vector<bool> used(n + 1, false);
        do {
            used[j0] = true;
            int i0 = p[j0];
            double delta = kInf;
            int j1 = 0;
            for (int j = 1; j <= n; ++j) {
                if (used[j]) continue;
                double c = std::abs(poles[j - 1] - prev[i0 - 1]) - u[i0] - v[j];
                if (c < minv[j]) { minv[j] = c; way[j] = j0; }
                if (minv[j] < delta) { delta = minv[j]; j1 = j; }
            }
            for (int j = 0; j <= n; ++j) {
                if (used[j]) { u[p[j]] += delta; v[j] -= delta; }
                else           minv[j] -= delta;
            }
            j0 = j1;
        } while (p[j0] != 0);
        do {
            int j1 = way[j0];
            p[j0] = p[j1];
            j0 = j1;
        } while (j0);
    }

    std::vector<int> perm(n);
    for (int j = 1; j <= n; ++j)
        if (p[j]) perm[p[j] - 1] = j - 1;
    return perm;
}
} // anonymous namespace

void matchPoles(std::vector<std::complex<double>>& poles,
                const std::vector<std::complex<double>>& prev) {
    int n = static_cast<int>(poles.size());
    std::vector<int> perm = hungarianAssign(poles, prev);
    std::vector<std::complex<double>> matched(n);
    for (int i = 0; i < n; ++i) matched[i] = poles[perm[i]];
    poles = matched;
}

PoleZeroResult computePoleZero(
    const LinearSystem& sys, int output_i, int input_j) {
    PoleZeroResult result;
    int n = sys.states();

    // Poles = eigenvalues of A.
    // n == 0 is a static gain: no poles, trivially stable.  Eigen::EigenSolver
    // on a 0x0 matrix is a real memory fault in BOTH build configs (SIGSEGV
    // under NDEBUG, SIGABRT in debug) — not a disabled sanity check — so this
    // guard is a correctness requirement.  See issue #5.
    result.is_stable = true;
    if (n > 0) {
        Eigen::EigenSolver<Eigen::MatrixXd> es(sys.A, false);
        result.poles.resize(n);
        for (int i = 0; i < n; ++i) {
            result.poles[i] = es.eigenvalues()(i);
            if (result.poles[i].real() >= -1e-10) {
                result.is_stable = false;
            }
        }
    }

    // Transmission zeros
    double D_ij = sys.D(output_i, input_j);

    if (std::abs(D_ij) > 1e-14) {
        // A static gain has no zeros; leave result.zeros empty rather than
        // pushing a 0x0 matrix into EigenSolver.  See issue #5.
        if (n > 0) {
            Eigen::MatrixXd A_z = sys.A -
                sys.B.col(input_j) * sys.C.row(output_i) / D_ij;
            Eigen::EigenSolver<Eigen::MatrixXd> zes(A_z, false);
            result.zeros.resize(n);
            for (int i = 0; i < n; ++i) {
                result.zeros[i] = zes.eigenvalues()(i);
            }
        }
    } else {
        Eigen::MatrixXd P = Eigen::MatrixXd::Zero(n + 1, n + 1);
        P.topLeftCorner(n, n) = sys.A;
        P.topRightCorner(n, 1) = sys.B.col(input_j);
        P.bottomLeftCorner(1, n) = sys.C.row(output_i);
        P(n, n) = D_ij;

        Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(n + 1, n + 1);
        Q.topLeftCorner(n, n) = Eigen::MatrixXd::Identity(n, n);

        Eigen::GeneralizedEigenSolver<Eigen::MatrixXd> ges(P, Q);
        auto alphas = ges.alphas();
        auto betas = ges.betas();

        for (int k = 0; k < n + 1; ++k) {
            if (std::abs(betas(k)) > 1e-10) {
                result.zeros.push_back(alphas(k) / betas(k));
            }
        }
    }

    return result;
}

std::vector<RootLocusPoint> computeRootLocus(
    const LinearSystem& sys, int output_i, int input_j,
    double k_min, double k_max, int num_points) {
    int n = sys.states();
    if (n == 0) return {};
    double D_ij = sys.D(output_i, input_j);
    std::vector<RootLocusPoint> result;
    result.reserve(num_points);

    for (int step = 0; step < num_points; ++step) {
        double K = k_min + (k_max - k_min) * step /
                   std::max(num_points - 1, 1);

        double denom = 1.0 + K * D_ij;
        double K_eff = (std::abs(denom) > 1e-12) ? K / denom : K;

        Eigen::MatrixXd A_cl =
            sys.A - K_eff * sys.B.col(input_j) * sys.C.row(output_i);
        Eigen::EigenSolver<Eigen::MatrixXd> es(A_cl, false);

        std::vector<std::complex<double>> poles(n);
        for (int i = 0; i < n; ++i) {
            poles[i] = es.eigenvalues()(i);
        }

        if (!result.empty()) {
            matchPoles(poles, result.back().poles);
        }

        result.push_back({K, poles});
    }

    return result;
}

std::vector<RootLocusPoint> computeStateFeedbackLocus(
    const LinearSystem& sys,
    const Eigen::MatrixXd& K,
    double alpha_min, double alpha_max, int num_points) {
    int n = sys.states();
    if (n == 0) return {};
    std::vector<RootLocusPoint> result;
    result.reserve(num_points);

    for (int step = 0; step < num_points; ++step) {
        double alpha = alpha_min + (alpha_max - alpha_min) * step /
                       std::max(num_points - 1, 1);

        Eigen::MatrixXd A_cl = sys.A - alpha * sys.B * K;
        Eigen::EigenSolver<Eigen::MatrixXd> es(A_cl, false);

        std::vector<std::complex<double>> poles(n);
        for (int i = 0; i < n; ++i) {
            poles[i] = es.eigenvalues()(i);
        }

        if (!result.empty()) {
            matchPoles(poles, result.back().poles);
        }

        result.push_back({alpha, poles});
    }

    return result;
}

}  // namespace caliburn
