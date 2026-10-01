// src/analysis/lqr.cpp
#include "lqr.h"
#include "system_properties.h"

#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <Eigen/LU>
#include <Eigen/SVD>
#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <utility>

namespace caliburn {
namespace {

// Sign-function iteration parameters.  The iteration converges quadratically
// and settles in 2-25 steps for every preset; 100 is a runaway guard, not a
// working budget.
constexpr double kSignTol = 1e-9;
constexpr int kMaxIter = 100;

// A closed-loop pole this close to the imaginary axis, relative to the size of
// the closed-loop matrix, is marginal rather than stable: the cost integral
// does not converge and calling the gain optimal would be a lie.  Rejecting it
// is deliberate, and is why the comparison below is not a bare `>= 0.0`.
constexpr double kStabTol = 1e-12;

// Q must be symmetric to this tolerance, and its smallest eigenvalue may dip
// this far below zero before the matrix counts as indefinite rather than as a
// positive-semidefinite matrix with rounding on its null space.
constexpr double kSymTol = 1e-12;
constexpr double kPsdTol = 1e-12;

bool isSymmetric(const Eigen::MatrixXd& M, double tol) {
    return M.rows() == M.cols()
        && (M - M.transpose()).norm() <= tol * std::max(1.0, M.norm());
}

// Solve A'P + PA - PBR^-1B'P + Q = 0 for the stabilizing P via the matrix sign
// function.  `R_chol` is the caller's already-validated Cholesky factor, so R
// is never factored twice.
//
// The method: the stable invariant subspace of the Hamiltonian carries P, and
// sign(H) is the projector that exposes it.  Newton's iteration for the sign
// function is Z <- (Z + Z^-1)/2, written below as Z - (Z - Z^-1)/2; the
// determinant scaling ck is the standard equilibration that stops the early
// iterates from wandering when H is badly scaled.
bool solveCARE(const Eigen::MatrixXd& A,
               const Eigen::MatrixXd& B,
               const Eigen::MatrixXd& Q,
               const Eigen::LLT<Eigen::MatrixXd>& R_chol,
               Eigen::MatrixXd& P) {
    const int n = static_cast<int>(A.rows());

    const Eigen::MatrixXd BRinvBt = B * R_chol.solve(B.transpose());
    Eigen::MatrixXd H(2 * n, 2 * n);
    H << A, BRinvBt,
         Q, -A.transpose();

    Eigen::MatrixXd Z = H;
    const double two_n = static_cast<double>(2 * n);
    bool converged = false;

    for (int iter = 0; iter < kMaxIter; ++iter) {
        const Eigen::MatrixXd Z_old = Z;

        // A singular iterate means the Hamiltonian has an eigenvalue on the
        // imaginary axis, so no stable/unstable splitting exists and the sign
        // function is undefined.  Bailing here beats dividing by zero and
        // returning a matrix of NaNs that later reads as a plausible gain.
        const double det = std::abs(Z.determinant());
        if (!std::isfinite(det) || det < std::numeric_limits<double>::min()) {
            return false;
        }

        Z *= std::pow(det, -1.0 / two_n);
        const Eigen::MatrixXd Z_inv = Z.inverse();
        Z = Z - 0.5 * (Z - Z_inv);
        if (!Z.allFinite()) return false;

        // Relative, matching the tolerances on Q and R above: an absolute
        // threshold reads as "not converged" on a badly scaled plant and stops
        // early on a finely scaled one.
        if ((Z - Z_old).norm() < kSignTol * std::max(1.0, Z.norm())) {
            converged = true;
            break;
        }
    }
    if (!converged) return false;

    // Z has converged to sign(H), whose n x n blocks are named W below.  P
    // solves [W12; W22 + I] P = [W11 + I; W21], an overdetermined 2n x n
    // system that is consistent in exact arithmetic;
    // the SVD gives the least-squares solution and tolerates rank deficiency
    // in the block, which a plain solve would not.
    const Eigen::MatrixXd W11 = Z.block(0, 0, n, n);
    const Eigen::MatrixXd W12 = Z.block(0, n, n, n);
    const Eigen::MatrixXd W21 = Z.block(n, 0, n, n);
    const Eigen::MatrixXd W22 = Z.block(n, n, n, n);

    const Eigen::MatrixXd I = Eigen::MatrixXd::Identity(n, n);
    Eigen::MatrixXd lhs(2 * n, n), rhs(2 * n, n);
    lhs << W12, W22 + I;
    rhs << W11 + I, W21;

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(
        lhs, Eigen::ComputeThinU | Eigen::ComputeThinV);
    P = svd.solve(rhs);
    if (!P.allFinite()) return false;

    // P is symmetric in exact arithmetic; the two triangles differ only by
    // rounding, and averaging them is what makes K = R^-1B'P exact rather than
    // dependent on which triangle it happened to read.
    P = 0.5 * (P + P.transpose()).eval();
    return true;
}

LqrResult failure(std::string why) {
    LqrResult r;
    r.error = std::move(why);
    return r;
}

}  // anonymous namespace

LqrResult computeLQR(const LinearSystem& sys,
                     const Eigen::MatrixXd& Q,
                     const Eigen::MatrixXd& R) {
    const int n = sys.states();
    const int m = sys.inputs();

    if (n == 0 || m == 0) return failure("plant has no states or no inputs");
    if (sys.A.cols() != n) return failure("A is not square");
    if (sys.B.rows() != n) return failure("B has the wrong number of rows");

    if (Q.rows() != n || Q.cols() != n) {
        return failure("Q must be " + std::to_string(n) + "x" +
                       std::to_string(n) + ", one entry per state");
    }
    if (R.rows() != m || R.cols() != m) {
        return failure("R must be " + std::to_string(m) + "x" +
                       std::to_string(m) + ", one entry per input");
    }

    if (!isSymmetric(Q, kSymTol)) return failure("Q must be symmetric");
    if (!isSymmetric(R, kSymTol)) return failure("R must be symmetric");

    // Q positive semidefinite: a negative state weight rewards deviation, and
    // the cost integral is then unbounded below.
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> q_eig(Q);
    if (q_eig.info() != Eigen::Success) {
        return failure("Q eigenvalue decomposition failed");
    }
    if (q_eig.eigenvalues().minCoeff() < -kPsdTol * std::max(1.0, Q.norm())) {
        return failure("Q must be positive semidefinite");
    }

    // R positive definite: R^-1 appears in the gain, so a singular R is an
    // input that costs nothing and would be driven without bound.  LLT is the
    // test as well as the factorization -- it fails exactly when R is not
    // positive definite.
    Eigen::LLT<Eigen::MatrixXd> R_chol(R);
    if (R_chol.info() != Eigen::Success) {
        return failure("R must be positive definite");
    }

    if (!isStabilizable(sys)) {
        return failure("plant is not stabilizable; an unstable mode is unreachable from the input");
    }

    Eigen::MatrixXd P;
    if (!solveCARE(sys.A, sys.B, Q, R_chol, P)) {
        return failure("Riccati solver did not converge");
    }

    LqrResult result;
    result.P = P;
    result.K = R_chol.solve(sys.B.transpose() * P);

    const Eigen::MatrixXd A_cl = sys.A - sys.B * result.K;
    Eigen::EigenSolver<Eigen::MatrixXd> es(A_cl);
    if (es.info() != Eigen::Success) {
        return failure("closed-loop eigenvalue decomposition failed");
    }
    result.closed_loop_poles.assign(
        es.eigenvalues().data(), es.eigenvalues().data() + n);

    // The stabilizing solution is the one the sign function is supposed to
    // pick.  If a closed-loop pole came back in the right half-plane it picked
    // the wrong invariant subspace, and the gain is worse than useless -- it
    // is a destabilizing gain wearing an optimal label.
    const double pole_scale = std::max(1.0, A_cl.norm());
    for (const std::complex<double>& pole : result.closed_loop_poles) {
        if (pole.real() > -kStabTol * pole_scale) {
            return failure("Riccati solution is not stabilizing");
        }
    }

    result.success = true;
    return result;
}

}  // namespace caliburn
