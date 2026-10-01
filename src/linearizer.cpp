// From caliburn reference/models/linearizer.cpp @ 08b44f2
#include "linearizer.h"

namespace caliburn {

LinearSystem linearize(
    const NonlinearFn& f,
    const Eigen::VectorXd& x0,
    const Eigen::VectorXd& u0,
    double epsilon)
{
    int n = static_cast<int>(x0.size());
    int m = static_cast<int>(u0.size());

    Eigen::MatrixXd C = Eigen::MatrixXd::Identity(n, n);
    Eigen::MatrixXd D = Eigen::MatrixXd::Zero(n, m);

    return linearize(f, x0, u0, C, D, epsilon);
}

LinearSystem linearize(
    const NonlinearFn& f,
    const Eigen::VectorXd& x0,
    const Eigen::VectorXd& u0,
    const Eigen::MatrixXd& C,
    const Eigen::MatrixXd& D,
    double epsilon)
{
    int n = static_cast<int>(x0.size());
    int m = static_cast<int>(u0.size());

    LinearSystem sys;
    sys.A = Eigen::MatrixXd::Zero(n, n);
    sys.B = Eigen::MatrixXd::Zero(n, m);
    sys.C = C;
    sys.D = D;

    // A = df/dx via central differences
    for (int i = 0; i < n; ++i) {
        Eigen::VectorXd x_plus = x0;
        Eigen::VectorXd x_minus = x0;
        x_plus(i) += epsilon;
        x_minus(i) -= epsilon;

        Eigen::VectorXd f_plus = f(x_plus, u0);
        Eigen::VectorXd f_minus = f(x_minus, u0);

        sys.A.col(i) = (f_plus - f_minus) / (2.0 * epsilon);
    }

    // B = df/du via central differences
    for (int j = 0; j < m; ++j) {
        Eigen::VectorXd u_plus = u0;
        Eigen::VectorXd u_minus = u0;
        u_plus(j) += epsilon;
        u_minus(j) -= epsilon;

        Eigen::VectorXd f_plus = f(x0, u_plus);
        Eigen::VectorXd f_minus = f(x0, u_minus);

        sys.B.col(j) = (f_plus - f_minus) / (2.0 * epsilon);
    }

    return sys;
}

ValidationResult validate(
    const LinearSystem& analytical,
    const NonlinearFn& f,
    const Eigen::VectorXd& x0,
    const Eigen::VectorXd& u0,
    double tol)
{
    auto numerical = linearize(f, x0, u0, analytical.C, analytical.D);

    ValidationResult result;
    result.A_error = (analytical.A - numerical.A).cwiseAbs();
    result.B_error = (analytical.B - numerical.B).cwiseAbs();
    result.max_A_error = result.A_error.maxCoeff();
    result.max_B_error = result.B_error.maxCoeff();
    result.pass = (result.max_A_error < tol) && (result.max_B_error < tol);

    return result;
}

}  // namespace caliburn
