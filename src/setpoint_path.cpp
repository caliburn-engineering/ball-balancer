// src/setpoint_path.cpp
#include "setpoint_path.h"

#include <algorithm>
#include <cmath>

namespace caliburn {
namespace {

/// Corners of the regular polygon a shape describes, or 0 for the others.
int corners(PathShape s) {
    switch (s) {
        case PathShape::Triangle: return 3;
        case PathShape::Square:   return 4;
        default:                  return 0;
    }
}

/// Corner `k` of an `n`-gon on the circumradius.
///
/// The first corner is put at the top rather than at +x so that the square
/// reads as a square: on this plate the visitor is looking at a triad of legs
/// at 120 degrees, and a square rotated 45 degrees against it looks like a
/// mistake rather than a choice.
Eigen::Vector2d corner(int k, int n, double r) {
    const double a = M_PI / 2.0 + 2.0 * M_PI * k / n;
    return Eigen::Vector2d(r * std::cos(a), r * std::sin(a));
}

/// Phase into [0, 1).  Callers are expected to hand one that already is, but a
/// phase is a lap position and a lap position wraps, so this is the definition
/// rather than a guard.
///
/// A non-finite phase cannot come out of `advancePhase`, but `pathPoint` is
/// public and the polygon walk below casts to `int` — which is undefined for a
/// NaN rather than merely wrong.  So it is stopped at the door.
double wrap01(double phase) {
    if (!std::isfinite(phase)) return 0.0;
    const double w = std::fmod(phase, 1.0);
    return (w < 0.0) ? w + 1.0 : w;
}

/// Position along the perimeter, as an edge and how far along it.
struct Along {
    int edge;
    double frac;   ///< 0..1 along that edge
};

Along along(double phase, int n) {
    const double travelled = wrap01(phase) * n;      // edges travelled
    const int edge = static_cast<int>(travelled) % n;
    return {edge, travelled - std::floor(travelled)};
}

}  // namespace

double advancePhase(double phase, double dt, double period_s) {
    if (period_s <= 0.0) return wrap01(phase);
    return wrap01(phase + dt / period_s);
}

Eigen::Vector2d pathPoint(const SetpointPath& p, double phase) {
    switch (p.shape) {
        case PathShape::Fixed:
            return Eigen::Vector2d::Zero();
        case PathShape::Circle: {
            const double a = 2.0 * M_PI * wrap01(phase);
            return Eigen::Vector2d(p.radius_m * std::cos(a), p.radius_m * std::sin(a));
        }
        default: break;
    }
    const int n = corners(p.shape);
    const Along a = along(phase, n);
    const Eigen::Vector2d from = corner(a.edge, n, p.radius_m);
    const Eigen::Vector2d to = corner(a.edge + 1, n, p.radius_m);
    return from + a.frac * (to - from);
}

Eigen::Vector2d pathVelocity(const SetpointPath& p, double phase) {
    if (p.shape == PathShape::Fixed || p.period_s <= 0.0)
        return Eigen::Vector2d::Zero();
    if (p.shape == PathShape::Circle) {
        const double a = 2.0 * M_PI * wrap01(phase);
        const double w = 2.0 * M_PI / p.period_s;     // rad/s
        return Eigen::Vector2d(-p.radius_m * w * std::sin(a),
                                p.radius_m * w * std::cos(a));
    }
    const int n = corners(p.shape);
    const Along a = along(phase, n);
    const Eigen::Vector2d from = corner(a.edge, n, p.radius_m);
    const Eigen::Vector2d to = corner(a.edge + 1, n, p.radius_m);
    // One edge per period/n, covered at constant speed.
    return (to - from) * (n / p.period_s);
}

double phaseNearest(const SetpointPath& p, const Eigen::Vector2d& target) {
    // The centre is equidistant from every point of any closed path round it,
    // so it has no nearest one and the answer is a choice rather than a
    // computation.  Made here, once, for every shape: without it the polygon
    // walk below picks whichever edge rounding happened to make shortest,
    // which is arbitrary AND liable to differ between builds.  Reachable in
    // one press — "Centre setpoint", then choose a shape.
    if (target.squaredNorm() <= 0.0) return 0.0;

    switch (p.shape) {
        case PathShape::Fixed:
            return 0.0;
        case PathShape::Circle:
            return wrap01(std::atan2(target(1), target(0)) / (2.0 * M_PI));
        default:
            break;
    }
    const int n = corners(p.shape);
    double best_d2 = -1.0, best_phase = 0.0;
    for (int k = 0; k < n; ++k) {
        const Eigen::Vector2d a = corner(k, n, p.radius_m);
        const Eigen::Vector2d edge = corner(k + 1, n, p.radius_m) - a;
        // Clamped, so the answer is a point ON the edge rather than on the
        // infinite line through it — which for a corner is the corner itself.
        const double len2 = edge.squaredNorm();
        const double f = (len2 > 0.0)
                             ? std::clamp((target - a).dot(edge) / len2, 0.0, 1.0)
                             : 0.0;
        const double d2 = (a + f * edge - target).squaredNorm();
        if (best_d2 < 0.0 || d2 < best_d2) {
            best_d2 = d2;
            best_phase = (k + f) / n;
        }
    }
    return wrap01(best_phase);
}

double clampPeriod(const SetpointPath& p, double asked_s) {
    return std::max(asked_s, minPeriod(p));
}

PathStep stepPath(const SetpointPath& p, double phase, double dt) {
    return {pathPoint(p, phase), pathVelocity(p, phase),
            advancePhase(phase, dt, p.period_s)};
}

double pathLength(const SetpointPath& p) {
    switch (p.shape) {
        case PathShape::Fixed:  return 0.0;
        case PathShape::Circle: return 2.0 * M_PI * p.radius_m;
        default: break;
    }
    const int n = corners(p.shape);
    return n * 2.0 * p.radius_m * std::sin(M_PI / n);
}

double minPeriod(const SetpointPath& p) {
    const double len = pathLength(p);
    if (len <= 0.0) return 0.0;
    return std::max(len / kMaxSetpointSpeed, kMinLapSeconds);
}

void pathOutline(const SetpointPath& p, int samples,
                 Eigen::Matrix<double, 2, Eigen::Dynamic>* out) {
    if (!out) return;
    if (p.shape == PathShape::Fixed) { out->resize(2, 0); return; }
    if (p.shape == PathShape::Circle) {
        out->resize(2, samples + 1);
        for (int i = 0; i <= samples; ++i) {
            const double a = 2.0 * M_PI * i / samples;
            out->col(i) = Eigen::Vector2d(p.radius_m * std::cos(a),
                                          p.radius_m * std::sin(a));
        }
        return;
    }
    const int n = corners(p.shape);
    out->resize(2, n + 1);
    for (int k = 0; k <= n; ++k) out->col(k) = corner(k, n, p.radius_m);
}

}  // namespace caliburn
