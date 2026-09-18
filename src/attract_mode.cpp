// src/attract_mode.cpp
#include "attract_mode.h"

namespace caliburn {

SetpointPath openingPath() {
    SetpointPath p;
    p.shape = PathShape::Circle;
    // radius_m and period_s take the type's own defaults — 120 mm and ten
    // seconds — so the opening and the sliders' starting position are one
    // fact rather than two that have to be kept in step.
    return p;
}

Eigen::Vector4d attractStart(const SetpointPath& path) {
    // Phase zero: the start of the lap, which is where `plate_view` opens its
    // accumulated phase.  The two have to agree or the demo's first frame is
    // already an error to answer.
    const Eigen::Vector2d p = pathPoint(path, 0.0);
    const Eigen::Vector2d v = pathVelocity(path, 0.0);
    return Eigen::Vector4d(p(0), p(1), v(0), v(1));
}

}  // namespace caliburn
