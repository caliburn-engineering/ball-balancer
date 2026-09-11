// src/auto_balance.cpp
#include "auto_balance.h"

#include <algorithm>
#include <cmath>

namespace caliburn {

namespace {
constexpr int kStates = 7;
constexpr int kLegs = 3;
}  // namespace

bool samePlant(const TableParams& a, double g_a,
               const TableParams& b, double g_b) {
    // A micron, not an epsilon.  The model's geometry arrives as `float`
    // sliders and the plate's is a `double` literal, so 0.300 and 0.300 differ
    // in the twelfth digit and an exact-ish comparison refuses two plates that
    // are the same object.  A micron is below anything the sliders can express
    // and far above the round-trip; the same argument sizes the gravity
    // tolerance, in m/s^2.
    const double tol = 1e-6;
    return std::abs(a.R_ground - b.R_ground) < tol &&
           std::abs(a.R_table - b.R_table) < tol &&
           std::abs(a.L1 - b.L1) < tol &&
           std::abs(a.L2 - b.L2) < tol &&
           std::abs(a.gamma_offset - b.gamma_offset) < tol &&
           std::abs(g_a - g_b) < tol;
}

bool gainFitsCascade(const AutoBalanceDesign& d) {
    return d.K.rows() == kLegs && d.K.cols() == kStates;
}

Eigen::VectorXd cascadeState(const std::array<double, 3>& alpha_rad,
                             double home_leg_rad,
                             const Eigen::Vector4d& ball) {
    Eigen::VectorXd x(kStates);
    for (int i = 0; i < kLegs; ++i) x(i) = alpha_rad[i] - home_leg_rad;
    x.segment<4>(kLegs) = ball;
    return x;
}

namespace {

/// `s` of the way from `from` toward `to`.
std::array<double, 3> lerp(const std::array<double, 3>& from,
                           const std::array<double, 3>& to, double s) {
    return {from[0] + s * (to[0] - from[0]),
            from[1] + s * (to[1] - from[1]),
            from[2] + s * (to[2] - from[2])};
}

}  // namespace

Retreat retreatToHoldable(const TableKinematics& tk,
                          const std::array<double, 3>& target,
                          const std::array<double, 3>& safe,
                          double condition_limit) {
    if (tk.can_hold(target, condition_limit)) return {target, false, false};
    // `safe` itself cannot be held: there is nothing to retreat TO.  Say so
    // rather than returning it as though it were a normal clip.
    if (!tk.can_hold(safe, condition_limit)) return {safe, true, true};

    // **March, then bisect inside the bracket** — the shape `max_conditioned_tilt`
    // already uses, and for the same reason.  A bisection over the whole ray
    // assumes the holdable points are one unbroken run from `safe`, and they
    // are not: holdability is a reachability set intersected with a condition
    // sublevel set, and the condition number rises toward a singularity and
    // falls again past it.  So a ray can run good, bad, good — measured over
    // 3000 random targets, 2 of the 1623 retreats bisection performed returned
    // a scale on the far side of an unholdable band, which is a command the
    // legs cannot travel to because the servo step is clipped at the band.
    //
    // The step is half a degree of the widest leg coordinate, so it costs what
    // the ray is long: one or two solves for a servo step, which is a fraction
    // of a degree and is what runs every frame, and up to seventy for a command
    // pinned to the far corner of the servo box, which is a frame that is
    // already saturated.  Half a degree is `max_conditioned_tilt`'s resolution
    // and is five times finer than the narrowest band that measurement found.
    constexpr double kMarchStep = 0.5 * M_PI / 180.0;
    constexpr int kRefinements = 12;   // 0.5 deg / 2^12, well past useful

    double span = 0.0;
    for (int i = 0; i < 3; ++i)
        span = std::max(span, std::abs(target[i] - safe[i]));
    const int steps = std::max(1, static_cast<int>(std::ceil(span / kMarchStep)));

    // `bad` starts at 1: `target` is already known unholdable, so the bracket
    // exists however far the march gets.
    double good = 0.0, bad = 1.0;
    for (int i = 1; i < steps; ++i) {
        const double s = static_cast<double>(i) / steps;
        if (tk.can_hold(lerp(safe, target, s), condition_limit)) good = s;
        else { bad = s; break; }
    }

    // Only now, inside a bracket the march has proved straddles an edge, is
    // halving safe.  `good` is only ever a scale that was TESTED holdable and
    // that the march reached by an unbroken run, so what comes back is both a
    // pose the plate can hold and one it can be driven to.
    for (int i = 0; i < kRefinements; ++i) {
        const double mid = 0.5 * (good + bad);
        if (tk.can_hold(lerp(safe, target, mid), condition_limit)) good = mid;
        else bad = mid;
    }
    return {lerp(safe, target, good), true, false};
}

LegCommand legCommand(const TableKinematics& tk,
                      const AutoBalanceDesign& d,
                      const std::array<double, 3>& alpha_rad,
                      const Eigen::Vector4d& ball,
                      const BallReference& ref) {
    LegCommand out{{d.home_leg_rad, d.home_leg_rad, d.home_leg_rad}, false, false};
    if (!gainFitsCascade(d)) return out;

    // The whole reference state: zero leg deviation, the setpoint, and the
    // speed the setpoint is travelling at.  The last two slots are the
    // velocity feedforward, and they are zero for a setpoint being held —
    // which is the case the loop spent its whole life in before #24.
    Eigen::VectorXd x_ref = Eigen::VectorXd::Zero(kStates);
    x_ref(3) = ref.position(0);
    x_ref(4) = ref.position(1);
    x_ref(5) = ref.velocity(0);
    x_ref(6) = ref.velocity(1);

    const Eigen::VectorXd u =
        -d.K * (cascadeState(alpha_rad, d.home_leg_rad, ball) - x_ref);

    for (int i = 0; i < kLegs; ++i) {
        const double raw = d.home_leg_rad + u(i);
        const double clamped = std::clamp(raw, d.alpha_min_rad, d.alpha_max_rad);
        if (clamped != raw) out.saturated = true;
        out.alpha_rad[i] = clamped;
    }

    // The travel limits are a box; the workspace is not, and the workspace's
    // own edge is not the bound either.  The level pose is the safe end here —
    // always assemblable AND far from any singularity, whatever the gain asked
    // for.
    const std::array<double, 3> level = {d.home_leg_rad, d.home_leg_rad,
                                         d.home_leg_rad};
    const Retreat r = retreatToHoldable(tk, out.alpha_rad, level,
                                        kRatesUntrustworthyAbove);
    out.alpha_rad = r.alpha_rad;
    out.clipped_to_holdable = r.retreated;
    return out;
}

std::array<double, 3> stepServos(const std::array<double, 3>& alpha_rad,
                                 const std::array<double, 3>& cmd_rad,
                                 double tau,
                                 double dt) {
    // tau <= 0 is a servo with no lag at all, which is what the plate had
    // before this existed.  It is reachable from the model panel's slider only
    // at its floor, but a degenerate tau must not divide.
    const double decay = (tau > 0.0) ? std::exp(-dt / tau) : 0.0;
    std::array<double, 3> next{};
    for (int i = 0; i < 3; ++i)
        next[i] = cmd_rad[i] + (alpha_rad[i] - cmd_rad[i]) * decay;
    return next;
}

Eigen::Vector2d predictedLanding(const Eigen::Matrix<double, 6, 1>& b,
                                 double ball_radius,
                                 double gravity) {
    const Eigen::Vector2d here(b(0), b(1));
    const double dz = b(2) - ball_radius;   // height above the surface
    if (dz <= 0.0 || gravity <= 0.0) return here;

    // z(t) = dz + vz t - g t^2 / 2 = 0, positive root.
    const double vz = b(5);
    const double disc = vz * vz + 2.0 * gravity * dz;
    if (disc < 0.0) return here;            // cannot happen for dz > 0, but say so
    const double t = (vz + std::sqrt(disc)) / gravity;
    return here + t * Eigen::Vector2d(b(3), b(4));
}

std::array<double, 3> stepServosOnPlate(const TableKinematics& tk,
                                        const std::array<double, 3>& alpha_rad,
                                        const std::array<double, 3>& cmd_rad,
                                        double tau,
                                        double dt) {
    // Where the legs ARE is the safe end: they started at home and have never
    // been moved anywhere the plate cannot be held, so it holds by induction.
    return retreatToHoldable(tk, stepServos(alpha_rad, cmd_rad, tau, dt),
                             alpha_rad, kRatesUntrustworthyAbove).alpha_rad;
}

namespace {

// The three tunings, measured rather than guessed.  Every number in the
// blurbs is asserted in `test_auto_balance`, so a claim on screen that stops
// being true fails a build rather than quietly misleading a visitor.
//
// R is 1 in all three.  Pricing the legs higher buys sluggishness the position
// weight already buys more legibly, and pricing them lower buys a tuning that
// loses the ball — measured, R = 0.1 against the nominal Q throws it off on 8
// of 72 kick directions.  One knob moves, so a visitor can say which one they
// just watched.
//
// The aggressive end is bounded by the ball, not by the solver.  Q on position
// at 300 recovers the 60/-40 mm displacement faster still, and then flings the
// ball clean off the plate on a 0.43 m/s nudge — the plate is moving out from
// under it rather than tilting under it (issue #23).  150 is the fastest
// tuning measured that loses the ball at no nudge magnitude the slider offers
// and the nominal tuning survives.
constexpr std::array<LqrPreset, 3> kPresets = {{
    {"Detuned",
     "sluggish - 11 s to settle against nominal's 1.8, and never near the stops",
     20.0, 2.0, 1.0},
    {"Nominal",
     "the tuning the demo opens on - 1.8 s, well damped, no saturation",
     100.0, 10.0, 1.0},
    {"Aggressive",
     "0.6 s, and living against the servo stops the whole way in",
     150.0, 2.0, 1.0},
}};

constexpr int kNominalIndex = 1;

}  // namespace

const std::array<LqrPreset, 3>& lqrPresets() { return kPresets; }

const LqrPreset& nominalPreset() { return kPresets[kNominalIndex]; }

Eigen::VectorXd presetStateWeights(const LqrPreset& p, int n) {
    Eigen::VectorXd q = Eigen::VectorXd::Ones(n);
    if (n == kStates) {
        q(3) = q(4) = p.q_position;  // ball position is what the product is about
        q(5) = q(6) = p.q_velocity;  // and enough velocity weight to arrive damped
    }
    return q;
}

Eigen::VectorXd presetInputWeights(const LqrPreset& p, int m) {
    return Eigen::VectorXd::Constant(m, p.r);
}

int activePreset(const Eigen::VectorXd& q, const Eigen::VectorXd& r) {
    // Off the cascade every preset collapses to unit weights, so a unit-weight
    // plant would "match" whichever one is listed first.  That is not a match,
    // it is the absence of a distinction, and saying -1 is the honest answer.
    if (q.size() != kStates || r.size() != kLegs) return -1;
    for (int i = 0; i < static_cast<int>(kPresets.size()); ++i) {
        if (q.isApprox(presetStateWeights(kPresets[i], kStates)) &&
            r.isApprox(presetInputWeights(kPresets[i], kLegs)))
            return i;
    }
    return -1;
}

Eigen::VectorXd defaultLqrStateWeights(int n) {
    return presetStateWeights(nominalPreset(), n);
}

Eigen::VectorXd defaultLqrInputWeights(int m) {
    return presetInputWeights(nominalPreset(), m);
}

}  // namespace caliburn
