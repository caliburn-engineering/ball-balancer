#pragma once

#include "table_kinematics.h"
#include <vector>
#include <functional>

/// Single point in a validated trajectory
struct TrajectoryPoint {
    double t;                       // Time [s]
    TablePose desired_pose;         // Commanded table pose
    IKResult ik;                    // Servo angles from IK
    std::array<double, 3> alpha_dot; // Servo angular velocities [rad/s]
    double condition_number;        // Velocity Jacobian condition number
    double manipulability;          // Manipulability index
    double det_J_pose;              // Determinant of J_pose (0 = singularity)

    bool pose_feasible;             // IK has valid solution within limits
    bool velocity_feasible;         // All servo velocities within limits
    bool well_conditioned;          // Condition number below threshold

    /// Overall: all constraints satisfied
    bool ok() const { return pose_feasible && velocity_feasible && well_conditioned; }
};

/// Validated trajectory result
struct TrajectoryResult {
    std::vector<TrajectoryPoint> points;
    int total_points;
    int feasible_count;
    int velocity_ok_count;
    int well_conditioned_count;
    double min_manipulability;
    double max_condition_number;
    double worst_det_J_pose;

    /// Index of first infeasible point (-1 if all feasible)
    int first_infeasible;
};

/// Configuration for trajectory validation
struct ValidatorConfig {
    /// Max servo angular velocity [rad/s].  The mechanism's own, not a second
    /// opinion about it: this was an uncited 10.0 and #32 gave the plant a
    /// cited 10.47, which would have left the validator rejecting leg rates the
    /// plate is happy to produce.  Nothing constructs a `ValidatorConfig` today
    /// — the first caller would have inherited the disagreement.
    double omega_max = kServoRateMax;
    double condition_threshold = 50.0; // Condition number warning threshold
    double dt = 0.001;                 // Timestep for finite differences [s]
};

class TrajectoryValidator {
public:
    TrajectoryValidator(const TableKinematics& tk, const ValidatorConfig& cfg = {});

    /// Validate a trajectory defined by a pose function over [t_start, t_end]
    TrajectoryResult validate(std::function<TablePose(double)> pose_fn,
                              double t_start, double t_end, double dt_sample) const;

    /// Validate a trajectory from a vector of (time, pose) pairs
    TrajectoryResult validate(const std::vector<std::pair<double, TablePose>>& trajectory) const;

    /// Print a summary of the validation result
    static void print_summary(const TrajectoryResult& result);

    /// Print detailed per-point data (for debugging / CSV export)
    static void print_detail(const TrajectoryResult& result);

private:
    const TableKinematics& tk_;
    ValidatorConfig cfg_;
};
