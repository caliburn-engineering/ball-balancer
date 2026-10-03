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
#include <string>
#include <utility>

namespace caliburn {
namespace {

// Sign-function iteration parameters.  The iteration converges quadratically
// and settles in 2-25 steps for every preset; 100 is a runaway guard, not a
// working budget.
constexpr double kSignTol = 1e-9;
constexpr int kMaxIter = 100;

// Kleinman refinement cap.  Quadratic convergence means 1-2 steps from a good
// sign-function starting point; 5 is insurance, not a working budget.
constexpr int kMaxRefine = 5;

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
// function is Z <- (Z + Z^-1)/2; norm scaling gamma = sqrt(||Z^-1||/||Z||)
// balances the two arms each step without computing a determinant.
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
    bool converged = false;

    for (int iter = 0; iter < kMaxIter; ++iter) {
        // Norm scaling: gamma = sqrt(||Z^-1|| / ||Z||) balances the two arms of
        // the Newton step without touching the determinant, which overflows or
        // underflows on a large-but-workable Hamiltonian.  If the inverse is
        // not finite the Hamiltonian has an eigenvalue on the imaginary axis and
        // the sign function is undefined.
        const Eigen::MatrixXd Z_inv = Z.inverse();
        if (!Z_inv.allFinite()) return false;
        const double zn = Z.norm();
        const double zi = Z_inv.norm();
        const double gamma = (zn > 0.0 && zi > 0.0) ? std::sqrt(zi / zn) : 1.0;
        const Eigen::MatrixXd Z_new = 0.5 * (gamma * Z + (1.0 / gamma) * Z_inv);
        if (!Z_new.allFinite()) return false;

        // Relative, matching the tolerances on Q and R above: an absolute
        // threshold reads as "not converged" on a badly scaled plant and stops
        // early on a finely scaled one.
        if ((Z_new - Z).norm() < kSignTol * std::max(1.0, Z_new.norm())) {
            Z = Z_new;
            converged = true;
            break;
        }
        Z = Z_new;
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

// Relative Riccati residual: equation norm over the sum of its terms' norms.
// Zero when all terms vanish.
double relRiccatiResidual(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B,
                          const Eigen::MatrixXd& Q,
                          const Eigen::LLT<Eigen::MatrixXd>& R_chol,
                          const Eigen::MatrixXd& P) {
    const Eigen::MatrixXd At_P    = A.transpose() * P;
    const Eigen::MatrixXd P_A     = P * A;
    const Eigen::MatrixXd PBRiBtP = P * B * R_chol.solve(B.transpose() * P);
    const Eigen::MatrixXd residual = At_P + P_A - PBRiBtP + Q;
    const double den = At_P.norm() + P_A.norm() + PBRiBtP.norm() + Q.norm();
    return (den > 0.0) ? residual.norm() / den : 0.0;
}

// Solve A'X + XA = C for symmetric X using the Kronecker product identity:
//   (I ⊗ A' + A^T ⊗ I) vec(X) = vec(C).
// The n^2 × n^2 system is cheap for the plant sizes this solver handles (n ≤ ~14).
// A must be stable for a unique positive-definite solution to exist when C is
// positive-definite.
bool solveLyapunov(const Eigen::MatrixXd& A, const Eigen::MatrixXd& C,
                   Eigen::MatrixXd& X) {
    const int n = static_cast<int>(A.rows());
    const int n2 = n * n;
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(n2, n2);
    const Eigen::MatrixXd At = A.transpose();

    // I ⊗ A': block-diagonal, each n×n block is A'.
    for (int i = 0; i < n; ++i)
        M.block(i * n, i * n, n, n) += At;

    // A^T ⊗ I: the (i,j)-th n×n block is A(j,i) * I_n.
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            M.block(i * n, j * n, n, n).diagonal().array() += A(j, i);

    const Eigen::VectorXd rhs = Eigen::Map<const Eigen::VectorXd>(C.data(), n2);
    const Eigen::FullPivLU<Eigen::MatrixXd> lu(M);
    if (!lu.isInvertible()) return false;
    const Eigen::VectorXd sol = lu.solve(rhs);
    if (!sol.allFinite()) return false;

    X = Eigen::Map<const Eigen::MatrixXd>(sol.data(), n, n);
    X = 0.5 * (X + X.transpose()).eval();
    return true;
}

// Newton/Kleinman refinement of a CARE solution.  Each step solves the
// Lyapunov equation with the current closed-loop A, giving quadratic
// convergence.  Stops when the relative residual stops falling or reaches
// near-machine-precision.  Returns the relative residual after refinement.
double kleinmanRefine(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B,
                      const Eigen::MatrixXd& Q,
                      const Eigen::LLT<Eigen::MatrixXd>& R_chol,
                      Eigen::MatrixXd& P) {
    double res = relRiccatiResidual(A, B, Q, R_chol, P);
    for (int iter = 0; iter < kMaxRefine; ++iter) {
        const Eigen::MatrixXd K    = R_chol.solve(B.transpose() * P);
        const Eigen::MatrixXd A_K  = A - B * K;
        // C = Q + K'RK = Q + P B R^{-1} B' P  (using symmetry of P and R)
        const Eigen::MatrixXd C    = Q + P * B * K;
        Eigen::MatrixXd P_new;
        if (!solveLyapunov(A_K, -C, P_new)) break;
        if (!P_new.allFinite()) break;
        const double new_res = relRiccatiResidual(A, B, Q, R_chol, P_new);
        if (new_res >= res) break;
        P   = P_new;
        res = new_res;
    }
    return res;
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

    // Controllability, per the design.  Stabilizability is the weaker condition
    // the CARE actually needs -- an uncontrollable but already-stable mode is
    // harmless -- but the check that exists here is the controllability one,
    // and rejecting a stabilizable-only plant is a false negative with a clear
    // message rather than a wrong gain.
    if (!checkControllability(sys).pass) {
        return failure("plant is not controllable; no stabilizing gain exists");
    }

    Eigen::MatrixXd P;
    if (!solveCARE(sys.A, sys.B, Q, R_chol, P)) {
        return failure("Riccati solver did not converge");
    }

    LqrResult result;
    result.pre_refinement_residual = relRiccatiResidual(sys.A, sys.B, Q, R_chol, P);
    result.residual = kleinmanRefine(sys.A, sys.B, Q, R_chol, P);
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
