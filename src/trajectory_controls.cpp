// src/trajectory_controls.cpp
#include "trajectory_controls.h"

namespace caliburn {

void TrajectoryControls::setShape(PathShape s, double& phase) {
    path_.shape = s;
    // The floor first, the nearest phase second — see the header.
    path_.period_s = clampPeriod(path_, path_.period_s);
    phase = phaseNearest(path_, setpoint_m_);
}

void TrajectoryControls::setSizeMm(double mm) {
    path_.radius_m = mm * 1e-3;
    path_.period_s = clampPeriod(path_, path_.period_s);
}

void TrajectoryControls::setLapS(double seconds) {
    path_.period_s = clampPeriod(path_, seconds);
}

void TrajectoryControls::setHeldSetpointMm(double x_mm, double y_mm) {
    setpoint_m_ = Eigen::Vector2d(x_mm * 1e-3, y_mm * 1e-3);
}

void TrajectoryControls::toInput(SimInput& in) const {
    in.path = path_;
    in.held_setpoint = setpoint_m_;
}

void TrajectoryControls::fromReport(const SimReport& r) {
    if (onAPath()) setpoint_m_ = r.setpoint;
}

}  // namespace caliburn
