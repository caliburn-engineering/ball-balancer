#include "ball_plant_linear.h"
#include "rolling_dynamics.h"
#include <memory>

namespace caliburn {

LinearSystem ballBalancerLinearModel(double gravity) {
    constexpr double k = 5.0 / 7.0;  // rolling factor, solid sphere

    LinearSystem sys;

    // State: [x, y, vx, vy], Input: [alpha, beta]
    sys.A = Eigen::MatrixXd::Zero(4, 4);
    sys.A(0, 2) = 1.0;  // dx/dt = vx
    sys.A(1, 3) = 1.0;  // dy/dt = vy

    sys.B = Eigen::MatrixXd::Zero(4, 2);
    // beta (Y-axis tilt) drives x-acceleration: sin(beta) ~ beta
    sys.B(2, 1) = k * gravity;
    // alpha (X-axis tilt) drives y-acceleration: sin(alpha) ~ alpha
    sys.B(3, 0) = k * gravity;

    // Output: measure x, y positions
    sys.C = Eigen::MatrixXd::Zero(2, 4);
    sys.C(0, 0) = 1.0;
    sys.C(1, 1) = 1.0;

    sys.D = Eigen::MatrixXd::Zero(2, 2);

    return sys;
}

NonlinearFn ballBalancerNonlinearFn(double ball_radius,
                                    double ball_mass,
                                    double rolling_friction,
                                    double plate_half_width,
                                    double gravity) {
    BallParams ball{ball_radius, ball_mass, rolling_friction};
    PlateParams plate{plate_half_width, gravity};
    auto dynamics = std::make_shared<RollingBallDynamics>(ball, plate);

    return [dynamics](const Eigen::VectorXd& state, const Eigen::VectorXd& input)
        -> Eigen::VectorXd {
        return dynamics->derivatives(
            Eigen::Vector4d(state),
            input(0),   // alpha
            input(1));  // beta
    };
}

}  // namespace caliburn
