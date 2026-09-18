#include "trajectory_validator.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <limits>

TrajectoryValidator::TrajectoryValidator(const TableKinematics& tk, const ValidatorConfig& cfg)
    : tk_(tk), cfg_(cfg) {}

TrajectoryResult TrajectoryValidator::validate(
    std::function<TablePose(double)> pose_fn,
    double t_start, double t_end, double dt_sample) const {

    std::vector<std::pair<double, TablePose>> trajectory;
    for (double t = t_start; t <= t_end + 1e-12; t += dt_sample) {
        trajectory.push_back({t, pose_fn(t)});
    }
    return validate(trajectory);
}

TrajectoryResult TrajectoryValidator::validate(
    const std::vector<std::pair<double, TablePose>>& trajectory) const {

    TrajectoryResult result;
    result.total_points = static_cast<int>(trajectory.size());
    result.feasible_count = 0;
    result.velocity_ok_count = 0;
    result.well_conditioned_count = 0;
    result.min_manipulability = std::numeric_limits<double>::infinity();
    result.max_condition_number = 0;
    result.worst_det_J_pose = std::numeric_limits<double>::infinity();
    result.first_infeasible = -1;

    result.points.reserve(trajectory.size());

    IKResult prev_ik{};
    bool has_prev = false;

    for (size_t idx = 0; idx < trajectory.size(); ++idx) {
        const auto& [t, desired] = trajectory[idx];
        TrajectoryPoint pt;
        pt.t = t;
        pt.desired_pose = desired;

        // --- IK feasibility ---
        pt.ik = tk_.inverse_kinematics(desired);
        pt.pose_feasible = pt.ik.feasible;

        if (pt.pose_feasible) {
            result.feasible_count++;

            // --- Servo velocities (finite difference) ---
            pt.velocity_feasible = true;
            if (has_prev) {
                double dt = t - trajectory[idx - 1].first;
                if (dt > 1e-12) {
                    for (int i = 0; i < 3; ++i) {
                        pt.alpha_dot[i] = (pt.ik.alpha[i] - prev_ik.alpha[i]) / dt;
                        if (std::abs(pt.alpha_dot[i]) > cfg_.omega_max) {
                            pt.velocity_feasible = false;
                        }
                    }
                } else {
                    pt.alpha_dot = {0, 0, 0};
                }
            } else {
                pt.alpha_dot = {0, 0, 0};
            }

            if (pt.velocity_feasible) result.velocity_ok_count++;

            // --- Jacobian analysis ---
            pt.condition_number = tk_.condition_number(pt.ik.alpha, desired);
            pt.manipulability = tk_.manipulability(pt.ik.alpha, desired);
            pt.det_J_pose = tk_.det_J_pose(pt.ik.alpha, desired);

            pt.well_conditioned = (pt.condition_number < cfg_.condition_threshold);
            if (pt.well_conditioned) result.well_conditioned_count++;

            // Track extremes
            result.min_manipulability = std::min(result.min_manipulability, pt.manipulability);
            result.max_condition_number = std::max(result.max_condition_number, pt.condition_number);
            result.worst_det_J_pose = std::min(result.worst_det_J_pose,
                                                std::abs(pt.det_J_pose));

            prev_ik = pt.ik;
            has_prev = true;
        } else {
            pt.alpha_dot = {0, 0, 0};
            pt.velocity_feasible = false;
            pt.condition_number = std::numeric_limits<double>::infinity();
            pt.manipulability = 0;
            pt.det_J_pose = 0;
            pt.well_conditioned = false;
            has_prev = false; // reset for velocity computation
        }

        if (!pt.ok() && result.first_infeasible == -1) {
            result.first_infeasible = static_cast<int>(idx);
        }

        result.points.push_back(pt);
    }

    return result;
}

void TrajectoryValidator::print_summary(const TrajectoryResult& result) {
    constexpr double DEG = M_PI / 180.0;

    std::cout << "=== Trajectory Validation Summary ===" << std::endl;
    std::cout << "  Total points:        " << result.total_points << std::endl;
    std::cout << "  Pose feasible:       " << result.feasible_count
              << "/" << result.total_points << std::endl;
    std::cout << "  Velocity feasible:   " << result.velocity_ok_count
              << "/" << result.total_points << std::endl;
    std::cout << "  Well-conditioned:    " << result.well_conditioned_count
              << "/" << result.total_points << std::endl;
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "  Min manipulability:  " << result.min_manipulability << std::endl;
    std::cout << "  Max condition #:     " << result.max_condition_number << std::endl;
    std::cout << "  Worst |det(J_pose)|: " << result.worst_det_J_pose << std::endl;

    if (result.first_infeasible >= 0) {
        const auto& pt = result.points[result.first_infeasible];
        std::cout << "  First issue at t=" << pt.t << "s: "
                  << "phi=" << pt.desired_pose.phi / DEG << "° "
                  << "theta=" << pt.desired_pose.theta / DEG << "°"
                  << std::endl;
        if (!pt.pose_feasible)
            std::cout << "    -> IK infeasible" << std::endl;
        if (!pt.velocity_feasible)
            std::cout << "    -> Servo velocity exceeded" << std::endl;
        if (!pt.well_conditioned)
            std::cout << "    -> Poorly conditioned (cond=" << pt.condition_number << ")" << std::endl;
    } else {
        std::cout << "  All points OK" << std::endl;
    }
}

void TrajectoryValidator::print_detail(const TrajectoryResult& result) {
    constexpr double DEG = M_PI / 180.0;

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "  t[s]   phi[°]  theta[°]  z[mm]   "
              << "a0[°]   a1[°]   a2[°]   "
              << "w0[°/s] w1[°/s] w2[°/s]  "
              << "cond    manip   status" << std::endl;
    std::cout << std::string(120, '-') << std::endl;

    for (const auto& pt : result.points) {
        std::cout << std::setw(6) << pt.t << "  "
                  << std::setw(6) << pt.desired_pose.phi / DEG << "  "
                  << std::setw(8) << pt.desired_pose.theta / DEG << "  "
                  << std::setw(6) << pt.desired_pose.z_c * 1000 << "  ";

        if (pt.pose_feasible) {
            std::cout << std::setw(6) << pt.ik.alpha[0] / DEG << "  "
                      << std::setw(6) << pt.ik.alpha[1] / DEG << "  "
                      << std::setw(6) << pt.ik.alpha[2] / DEG << "  "
                      << std::setw(7) << pt.alpha_dot[0] / DEG << " "
                      << std::setw(7) << pt.alpha_dot[1] / DEG << " "
                      << std::setw(7) << pt.alpha_dot[2] / DEG << "  "
                      << std::setw(7) << pt.condition_number << " "
                      << std::setw(7) << pt.manipulability << "  ";
        } else {
            std::cout << "  ---     ---     ---     "
                      << "  ---     ---     ---     "
                      << "  ---     ---    ";
        }

        if (pt.ok()) std::cout << "OK";
        else {
            if (!pt.pose_feasible) std::cout << "INFEAS ";
            if (!pt.velocity_feasible) std::cout << "VIOL ";
            if (!pt.well_conditioned) std::cout << "COND ";
        }
        std::cout << std::endl;
    }
}
