// src/setpoint_path.cpp
#include "setpoint_path.h"

#include <algorithm>
#include <cmath>

namespace caliburn {
namespace {

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

/// A polygon with its corners filleted, measured out once so that every
/// question about it is asked of the same arithmetic.
///
/// The path is `n` identical units laid end to end, and a unit is half a
/// fillet, a straight, and half the next fillet.  Phase `k/n` sits at the
/// middle of corner `k`'s fillet, where the sharp polygon's corner `k` was —
/// so the phase means the same place it always did, and `phaseNearest`, the
/// size slider and the shape combo all keep the promises they made.
struct Filleted {
    int n = 0;              ///< corners; zero for the shapes that have none
    double R = 0.0;         ///< circumradius of the sharp polygon
    double half_angle = 0.0;///< pi/n: half the turn at each corner
    double side = 0.0;      ///< the sharp side length
    double rho = 0.0;       ///< fillet radius
    double inset = 0.0;     ///< centre to fillet centre: R - rho/cos(pi/n)
    double half_arc = 0.0;  ///< half of one fillet's arc length
    double straight = 0.0;  ///< what is left of a side between two fillets
    double unit = 0.0;      ///< perimeter per corner
    double perimeter = 0.0;

    /// The direction from the centre to corner `k`, as an angle.  The first
    /// corner is at the top — see `corner`.
    double corner_angle(int k) const {
        return M_PI / 2.0 + 2.0 * half_angle * k;
    }
    /// The centre of corner `k`'s fillet.
    Eigen::Vector2d fillet_centre(int k) const {
        const double a = corner_angle(k);
        return Eigen::Vector2d(inset * std::cos(a), inset * std::sin(a));
    }
    /// A point on corner `k`'s fillet circle at angle `chi`, and the unit
    /// tangent there.  Anticlockwise, which is the way the lap runs.
    Eigen::Vector2d on_fillet(int k, double chi) const {
        return fillet_centre(k) + rho * Eigen::Vector2d(std::cos(chi), std::sin(chi));
    }
    static Eigen::Vector2d tangent_at(double chi) {
        return Eigen::Vector2d(-std::sin(chi), std::cos(chi));
    }
};

/// The side of a regular `n`-gon on circumradius `R`, and the perimeter it has
/// once a fillet of radius `rho` is cut into each of its corners.
///
/// `n*side - rho*(2n*tan(pi/n) - 2*pi)`: each corner gives up two tangent
/// lengths and gets back a shorter arc.  Here rather than at each of the three
/// places that want it — the fillet solve, the geometry, and `minPeriod`,
/// which cannot use the other two because it is solving for the lap they take
/// as given.
double polygonSide(int n, double R) {
    return 2.0 * R * std::sin(M_PI / n);
}

double filletedPerimeter(int n, double side, double rho) {
    return n * side - rho * (2.0 * n * std::tan(M_PI / n) - 2.0 * M_PI);
}

/// The fillet radius a polygon's corners take at the speed it is being walked.
///
/// `rho = v^2 / a_max` with `v = perimeter / period`, and the perimeter is
/// itself `n*side - B*rho` — so this is circular, and it is solved rather than
/// iterated.  Writing `A = n*side` and `B = 2n*tan(pi/n) - 2*pi`:
///
///     rho * T^2 * a_max = (A - B*rho)^2
///     B^2 rho^2 - (2AB + T^2 a_max) rho + A^2 = 0
///
/// with the smaller root the one that is a fillet — the larger one is the
/// arithmetic's other branch, a "fillet" bigger than the shape.  Its
/// discriminant is `c*(c + 4AB)` with `c = T^2 a_max`, which is never negative
/// for a real path, so there is no case where the fillet fails to exist.
///
/// `B > 0` for every `n >= 3`: it is 4.11 for a triangle and 1.72 for a
/// square, and it is the length a unit fillet SAVES — the two tangent lengths
/// it removes are longer than the arc it puts back.
double filletFor(int n, double side, double period_s, double accel_max) {
    const double A = n * side;
    const double B = A - filletedPerimeter(n, side, 1.0);   // the per-metre saving
    const double c = period_s * period_s * accel_max;
    if (!(c > 0.0) || !(A > 0.0) || !(B > 0.0) || !std::isfinite(c)) return 0.0;
    const double b = 2.0 * A * B + c;
    const double rho = (b - std::sqrt(c * (c + 4.0 * A * B))) / (2.0 * B * B);
    return std::isfinite(rho) ? std::max(rho, 0.0) : 0.0;
}

/// Everything about the filleted polygon `p` describes.  Cheap enough to build
/// per call, and building it per call is what keeps the size slider, the lap
/// slider and the plant's `accel_max` from ever being read at different
/// instants.
Filleted filletOf(const SetpointPath& p) {
    Filleted f;
    f.n = pathCorners(p.shape);
    if (f.n == 0) return f;
    f.R = std::max(p.radius_m, 0.0);
    f.half_angle = M_PI / f.n;
    f.side = polygonSide(f.n, f.R);

    // Capped at the inradius, where the fillets meet each other and the polygon
    // has become its own incircle.  Past that a fillet would have to eat into
    // the previous corner, which is not a rounded polygon, it is a smaller one.
    const double inradius = f.R * std::cos(f.half_angle);
    f.rho = std::min(filletFor(f.n, f.side, p.period_s, p.accel_max), inradius);

    f.inset = f.R - ((f.rho > 0.0) ? f.rho / std::cos(f.half_angle) : 0.0);
    f.half_arc = f.rho * f.half_angle;
    f.straight = std::max(f.side - 2.0 * f.rho * std::tan(f.half_angle), 0.0);
    f.unit = 2.0 * f.half_arc + f.straight;
    f.perimeter = f.n * f.unit;
    return f;
}

/// Where the setpoint is at `phase`, and which way it is going.
///
/// The unit is walked in three pieces and the tangent falls out of each of
/// them, which is the point of filleting the GEOMETRY rather than slewing the
/// velocity: there is one curve here and this is a point on it and its
/// derivative, so `v = dp/dt` holds by construction instead of by agreement.
///
/// With `rho == 0` the two arc pieces have zero length and are never entered,
/// the fillet centre is the corner itself, and this is exactly the sharp walk
/// it replaced — same phase convention, same points.
struct Riding {
    Eigen::Vector2d point{Eigen::Vector2d::Zero()};
    Eigen::Vector2d tangent{Eigen::Vector2d::Zero()};   ///< unit, or zero
};

Riding ride(const Filleted& f, double phase) {
    Riding r;
    if (f.n == 0) return r;
    if (!(f.unit > 0.0)) {                      // a path of no size
        r.point = corner(0, f.n, f.R);
        return r;
    }
    const double travelled = wrap01(phase) * f.n;
    const int k = static_cast<int>(travelled) % f.n;
    const double local = (travelled - std::floor(travelled)) * f.unit;

    if (local < f.half_arc) {
        // Leaving corner k's fillet: from its middle towards the next side.
        const double chi = f.corner_angle(k) + local / f.rho;
        r.point = f.on_fillet(k, chi);
        r.tangent = Filleted::tangent_at(chi);
        return r;
    }
    // The fillet's outgoing tangent point, and the direction the side runs in.
    // Read off the fillet circle rather than off the corners, so that the two
    // agree to the last bit where they meet.
    const double chi_out = f.corner_angle(k) + f.half_angle;
    if (local < f.half_arc + f.straight) {
        r.tangent = Filleted::tangent_at(chi_out);
        r.point = f.on_fillet(k, chi_out) + (local - f.half_arc) * r.tangent;
        return r;
    }
    // Arriving at corner k+1's fillet, up to its middle.
    const double chi = chi_out + (local - f.half_arc - f.straight) / f.rho;
    r.point = f.on_fillet(k + 1, chi);
    r.tangent = Filleted::tangent_at(chi);
    return r;
}

}  // namespace

int pathCorners(PathShape s) {
    switch (s) {
        case PathShape::Triangle: return 3;
        case PathShape::Square:   return 4;
        default:                  return 0;
    }
}

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
    return ride(filletOf(p), phase).point;
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
    // Constant speed all the way round, fillets included: the whole perimeter
    // in one period.  The direction is the curve's own tangent, so this is the
    // position's time derivative rather than a second opinion about it.
    const Filleted f = filletOf(p);
    return ride(f, phase).tangent * (f.perimeter / p.period_s);
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
    const Filleted f = filletOf(p);
    if (!(f.unit > 0.0)) return 0.0;

    // Distance along the lap from phase zero — the middle of corner 0's fillet
    // — for each candidate, converted to a phase at the end.  The two kinds of
    // piece are measured the same way so that the comparison between them is
    // one comparison rather than two conventions.
    double best_d2 = -1.0, best_s = 0.0;
    const auto offer = [&](const Eigen::Vector2d& at, double s) {
        const double d2 = (at - target).squaredNorm();
        if (best_d2 < 0.0 || d2 < best_d2) { best_d2 = d2; best_s = s; }
    };

    for (int k = 0; k < f.n; ++k) {
        // The straight leaving corner k, clamped so the answer is a point ON
        // it rather than on the infinite line through it.
        const double chi_out = f.corner_angle(k) + f.half_angle;
        const Eigen::Vector2d from = f.on_fillet(k, chi_out);
        const Eigen::Vector2d dir = Filleted::tangent_at(chi_out);
        const double d = std::clamp((target - from).dot(dir), 0.0, f.straight);
        offer(from + d * dir, k * f.unit + f.half_arc + d);

        // And corner k's fillet, clamped in ANGLE about the fillet's centre —
        // the same projection one turn round, and the reason a fillet cannot be
        // approximated by its chord here: the nearest point on an arc is very
        // often in the middle of it.
        if (f.rho > 0.0) {
            const Eigen::Vector2d from_centre = target - f.fillet_centre(k);
            const double mid = f.corner_angle(k);
            // Into (-pi, pi] about the fillet's middle before clamping, so that
            // a target on the far side of the plate is behind the fillet rather
            // than a full turn ahead of it.
            double delta = 0.0;
            if (from_centre.squaredNorm() > 0.0) {
                delta = std::atan2(from_centre(1), from_centre(0)) - mid;
                delta = std::remainder(delta, 2.0 * M_PI);
            }
            delta = std::clamp(delta, -f.half_angle, f.half_angle);
            offer(f.on_fillet(k, mid + delta), k * f.unit + delta * f.rho);
        }
    }
    return wrap01(best_s / f.perimeter);
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
    return filletOf(p).perimeter;
}

double filletRadius(const SetpointPath& p) {
    return filletOf(p).rho;
}

double minPeriod(const SetpointPath& p) {
    const int n = pathCorners(p.shape);
    if (n == 0) {
        const double len = pathLength(p);
        if (len <= 0.0) return 0.0;
        return std::max(len / kMaxSetpointSpeed, kMinLapSeconds);
    }
    // A polygon's length depends on its lap, so `pathLength` would answer
    // about the lap the caller happens to hold rather than the one being
    // solved for.  At the cap the speed is known, though, and that breaks the
    // circle: the fillet is `kMaxSetpointSpeed^2 / accel_max` whatever lap
    // comes out, so the perimeter is, so the lap that walks that perimeter at
    // the cap is.  Substituting it back reproduces itself — it is the fixed
    // point, found in closed form.
    const double R = std::max(p.radius_m, 0.0);
    const double side = polygonSide(n, R);
    if (side <= 0.0) return 0.0;

    const double rho = (p.accel_max > 0.0)
                           ? std::min(kMaxSetpointSpeed * kMaxSetpointSpeed / p.accel_max,
                                      R * std::cos(M_PI / n))   // the inradius cap
                           : 0.0;
    return std::max(filletedPerimeter(n, side, rho) / kMaxSetpointSpeed,
                    kMinLapSeconds);
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
    const Filleted f = filletOf(p);
    if (!(f.rho > 0.0)) {                       // the sharp polygon, as before
        out->resize(2, f.n + 1);
        for (int k = 0; k <= f.n; ++k) out->col(k) = corner(k, f.n, f.R);
        return;
    }
    // Every fillet gets the same share of `samples`, and the straights need no
    // samples at all — a segment between two arcs' endpoints IS the straight.
    // At least two segments per fillet, so that a fillet is never drawn as the
    // chord it is there to replace.
    const int per_fillet = std::max(2, samples / (2 * f.n));
    out->resize(2, f.n * (per_fillet + 1) + 1);
    int col = 0;
    for (int k = 0; k < f.n; ++k) {
        const double mid = f.corner_angle(k);
        for (int i = 0; i <= per_fillet; ++i) {
            const double chi =
                mid + f.half_angle * (2.0 * i / per_fillet - 1.0);
            out->col(col++) = f.on_fillet(k, chi);
        }
    }
    out->col(col) = out->col(0);                // closed
}

}  // namespace caliburn
