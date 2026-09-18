// src/ball_contact.cpp
#include "ball_contact.h"

#include <algorithm>
#include <cmath>

namespace caliburn {
namespace {

/// The axial vector of a skew-symmetric matrix.
Eigen::Vector3d vee(const Eigen::Matrix3d& S) {
    return Eigen::Vector3d(S(2, 1), S(0, 2), S(1, 0));
}

}  // namespace

std::array<double, 3> servoRate(const std::array<double, 3>& alpha_rad,
                                const std::array<double, 3>& cmd_rad,
                                double tau,
                                double rate_max) {
    const bool limited = rateLimited(rate_max);
    std::array<double, 3> out{};
    for (int i = 0; i < 3; ++i) {
        const double err = cmd_rad[i] - alpha_rad[i];
        if (tau > 0.0) {
            const double lag = err / tau;
            out[i] = limited ? std::clamp(lag, -rate_max, rate_max) : lag;
        } else if (limited && err != 0.0) {
            // No lag to slow it, so the limit is the whole model.
            out[i] = (err > 0.0) ? rate_max : -rate_max;
        }
    }
    return out;
}

std::array<double, 3> servoAccel(const std::array<double, 3>& alpha_dot_rad_s,
                                 double tau,
                                 double rate_max) {
    const bool limited = rateLimited(rate_max);
    std::array<double, 3> out{};
    for (int i = 0; i < 3; ++i) {
        const double rate = alpha_dot_rad_s[i];
        // `>=`, not `>`, and the difference is load-bearing: `servoRate`
        // CLAMPS, so a leg that is genuinely ramping reports exactly
        // `rate_max` and nothing above it.  `>` would find no saturated leg
        // anywhere and quietly delete the property this exists for.
        if (limited && std::abs(rate) >= rate_max) continue;   // ramping: zero
        if (tau <= 0.0) continue;
        out[i] = -rate / tau;
    }
    return out;
}

PlateMotion plateMotion(const TableKinematics& tk,
                        const TablePose& pose,
                        const std::array<double, 3>& alpha_rad,
                        const std::array<double, 3>& alpha_dot_rad_s,
                        const std::array<double, 3>& alpha_ddot_rad_s2,
                        const PlateMotion* prev,
                        double dt) {
    PlateMotion m;
    m.R = tk.table_rotation(pose.phi, pose.theta);
    m.c = Eigen::Vector3d(0.0, 0.0, pose.z_c);

    // d(pose)/dt = J_v * d(alpha)/dt, with pose = (phi, theta, z_c).
    const Eigen::Vector3d alpha_dot(alpha_dot_rad_s[0], alpha_dot_rad_s[1],
                                    alpha_dot_rad_s[2]);
    const Eigen::Vector3d alpha_ddot(alpha_ddot_rad_s2[0], alpha_ddot_rad_s2[1],
                                     alpha_ddot_rad_s2[2]);
    m.J_v = tk.velocity_jacobian(alpha_rad, pose);
    const Eigen::Vector3d pose_dot = m.J_v * alpha_dot;
    m.c_dot = Eigen::Vector3d(0.0, 0.0, pose_dot(2));

    // omega is LINEAR in the tilt rates, so it factors as `A(pose) * (phi_dot,
    // theta_dot)` and the map can be read off by evaluating it at the two unit
    // rates.  Extracting it costs two rotation derivatives and buys the second
    // derivative below: differentiating a product of a smooth matrix and a
    // discontinuous vector is only possible once they are separated.
    m.A.col(0) = vee(tk.table_rotation_dot(pose.phi, pose.theta, 1.0, 0.0) *
                     m.R.transpose());
    m.A.col(1) = vee(tk.table_rotation_dot(pose.phi, pose.theta, 0.0, 1.0) *
                     m.R.transpose());
    m.omega = m.A * pose_dot.head<2>();

    // pose_ddot = J_v alpha_ddot + J_v_dot alpha_dot, and omega_dot the same
    // shape one level up.  The analytic halves carry the command step; the
    // differenced halves carry only how the geometry itself is changing, which
    // is continuous because the leg angles and the pose are.
    Eigen::Vector3d pose_ddot = m.J_v * alpha_ddot;
    if (prev != nullptr && dt > 0.0)
        pose_ddot += ((m.J_v - prev->J_v) / dt) * alpha_dot;
    m.c_ddot = Eigen::Vector3d(0.0, 0.0, pose_ddot(2));

    m.omega_dot = m.A * pose_ddot.head<2>();
    if (prev != nullptr && dt > 0.0)
        m.omega_dot += ((m.A - prev->A) / dt) * pose_dot.head<2>();

    m.rates_trustworthy =
        tk.condition_number(alpha_rad, pose) < kRatesUntrustworthyAbove;
    return m;
}

double normalAccel(const PlateMotion& plate,
                   const Eigen::Vector3d& s,
                   const Eigen::Vector2d& s_dot,
                   double gravity) {
    const Eigen::Vector3d n = plate.normal();
    const Eigen::Vector3d r_world = plate.R * s;            // centre, from the table centre
    const Eigen::Vector3d v_world = plate.R * Eigen::Vector3d(s_dot(0), s_dot(1), 0.0);

    const double gravity_share = gravity * n.z();
    const double heave = n.dot(plate.c_ddot);
    const double angular = n.dot(plate.omega_dot.cross(r_world) +
                                 plate.omega.cross(plate.omega.cross(r_world)));
    const double coriolis = 2.0 * n.dot(plate.omega.cross(v_world));

    return gravity_share + heave + angular + coriolis;
}

double quasiStaticNormalAccel(const PlateMotion& plate, double gravity) {
    return gravity * plate.normal().z();
}

double contactNormalRate(const PlateMotion& plate, const Eigen::Vector3d& s) {
    const Eigen::Vector3d n = plate.normal();
    return n.dot(plate.c_dot + plate.omega.cross(plate.R * s));
}

void worldOf(const BallState& b, const PlateMotion& plate, double ball_radius,
             Eigen::Vector3d* p, Eigen::Vector3d* v) {
    if (b.airborne) {
        if (p) *p = b.flight_p;
        if (v) *v = b.flight_v;
        return;
    }
    const Eigen::Vector3d s(b.rolling(0), b.rolling(1), ball_radius);
    const Eigen::Vector3d r_world = plate.R * s;
    if (p) *p = plate.c + r_world;
    // The three parts of a carried ball's velocity: the table centre moving,
    // the table rotating about it, and the ball rolling across it.  Drop the
    // first two and a ball launched off a slamming plate leaves with only its
    // rolling speed, which is the smallest part of the answer.
    if (v) {
        *v = plate.c_dot + plate.omega.cross(r_world) +
             plate.R * Eigen::Vector3d(b.rolling(2), b.rolling(3), 0.0);
    }
}

Eigen::Matrix<double, 6, 1> plateFrame(const BallState& b,
                                       const PlateMotion& plate,
                                       double ball_radius) {
    Eigen::Matrix<double, 6, 1> out;
    if (!b.airborne) {
        out << b.rolling(0), b.rolling(1), ball_radius,
               b.rolling(2), b.rolling(3), 0.0;
        return out;
    }
    // Position is a plain change of basis.  Velocity is not: the plate frame
    // is moving, so what it sees is the ball's world velocity minus the
    // velocity of the plate point the ball is currently above.
    const Eigen::Vector3d rel = b.flight_p - plate.c;
    const Eigen::Vector3d q = plate.R.transpose() * rel;
    const Eigen::Vector3d v_rel =
        b.flight_v - plate.c_dot - plate.omega.cross(rel);
    const Eigen::Vector3d qd = plate.R.transpose() * v_rel;
    out << q(0), q(1), q(2), qd(0), qd(1), qd(2);
    return out;
}

BallState stepBallContact(const RollingBallDynamics& dynamics,
                          const BallState& b,
                          const PlateMotion& now,
                          const PlateMotion& prev,
                          const TablePose& pose,
                          double ball_radius,
                          double gravity,
                          double dt,
                          ContactReport* report) {
    BallState out = b;
    // Zeroed first: a ball already in flight is pressed by nothing, and "no
    // impact this frame" is every frame but the arrivals.  Both are overwritten
    // by whichever branch below applies.
    if (report) *report = ContactReport{};

    if (!out.airborne) {
        // No opinion without trustworthy rates.  Staying in contact is the
        // conservative answer and the honest one: the ball was on the plate a
        // moment ago and nothing believable says it left.
        //
        // BOTH frames have to be trustworthy, not just this one.  `now`'s
        // accelerations are analytic in the servo lag but still carry a
        // differenced term — how `J_v` and `A` are themselves changing — and
        // `J_v` is `-J_pose^-1 J_alpha`, the very thing that blows up near a
        // singularity.  So a believable frame differenced against a garbage one
        // is still a garbage acceleration, and guarding only `now` would let
        // the first good frame on the way out of a singularity report a
        // separation built entirely from the rates the guard exists to refuse.
        if (!now.rates_trustworthy || !prev.rates_trustworthy) {
            const double N_qs = quasiStaticNormalAccel(now, gravity);
            if (report) report->normal_accel = N_qs;
            out.rolling = stepBall(dynamics, out.rolling, pose, N_qs, dt);
            return out;
        }
        const double N = normalAccel(now,
                                     Eigen::Vector3d(out.rolling(0), out.rolling(1),
                                                     ball_radius),
                                     Eigen::Vector2d(out.rolling(2), out.rolling(3)),
                                     gravity);
        if (report) report->normal_accel = N;
        if (N > 0.0) {
            // The same N that decided the ball is still touching also says how
            // hard it is pressed, and rolling resistance is a normal-force
            // effect: a ball on a plate dropping away beneath it is barely
            // retarded at all, and one in a plate's upswing is retarded more
            // than its weight alone would explain.  #23.
            out.rolling = stepBall(dynamics, out.rolling, pose, N, dt);
            return out;
        }
        // Contact is over.  The ball keeps the velocity it had, all of it.
        //
        // Converted against `prev`, not `now`: the ball was on the plate at the
        // START of this frame, so that is where it was and how fast it was
        // going.  Using the end-of-frame plate teleports the ball down by
        // however far the plate fell during the frame, and then the landing
        // test below finds it already touching — which reads as a ball that
        // cannot leave at all.
        worldOf(out, prev, ball_radius, &out.flight_p, &out.flight_v);
        out.airborne = true;
    }

    // Ballistic: gravity and nothing else, integrated exactly because a
    // constant acceleration has a closed form and there is no reason to ask
    // RK4 for an answer that is already written down.
    const Eigen::Vector3d g(0.0, 0.0, -gravity);
    const BallState open = out;   // the frame as it opened, for the rewind below
    const double gap_open =
        plateFrame(open, prev, ball_radius)(2) - ball_radius;

    out.flight_p += out.flight_v * dt + 0.5 * g * dt * dt;
    out.flight_v += g * dt;

    // Arrived?  Against `now`, the end-of-frame plate — the ball has had its
    // whole step and so has the plate, so this asks where they both ended up.
    // Measured in the plate's frame, so a plate that has tilted or heaved up to
    // meet the ball counts as catching it.
    const double gap_close =
        plateFrame(out, now, ball_radius)(2) - ball_radius;
    if (gap_close > 0.0) return out;

    // **The impact is somewhere inside this frame, and where exactly matters.**
    //
    // Bouncing at the frame boundary instead does not merely blur the answer,
    // it reverses its sign: the ball is found already BELOW the surface, having
    // fallen past it since the last sample, so the approach speed read there is
    // larger than the speed it truly arrived at — by up to `g dt`.  Reflecting
    // that at `e` gives back more than was brought in, and the sampling pumps
    // the bounce train instead of damping it.  Measured before this was fixed,
    // a train started at 1 m/s climbed to a fixed point near 1.3 m/s and never
    // terminated at all.
    //
    // So the crossing is found first.  The gap closes because the ball falls
    // and because the plate rises, and both are already in these two numbers:
    // `gap_open` is measured against the plate at the frame's start and
    // `gap_close` against the plate at its end.  Interpolating linearly between
    // them errs on the safe side — the true gap is concave under gravity, so
    // the chord crosses zero no later than the curve does, and the bounce is
    // taken at or before the real contact, never after it.  The residual is a
    // slight loss rather than a gain, which is the direction a discretisation
    // has to err in if the train is to end.
    //
    // A ball flush with the plate — `gap_open == 0`, which is exactly how a
    // ball that has only just separated opens its first airborne frame — makes
    // the interpolation degenerate, and the sign of the relative normal
    // velocity is what settles it instead.  A ball that is not moving into the
    // surface at the frame's open crosses later in the frame, and the end is
    // the earliest instant there is evidence for; a ball already moving into it
    // is in contact at the open.
    const double crossing =
        (gap_open > 0.0)
            ? gap_open / (gap_open - gap_close)
            : (plateFrame(open, prev, ball_radius)(5) < 0.0 ? 0.0 : 1.0);
    const double t_hit = std::clamp(crossing, 0.0, 1.0) * dt;

    BallState hit = open;
    hit.flight_p += open.flight_v * t_hit + 0.5 * g * t_hit * t_hit;
    hit.flight_v += g * t_hit;

    // Read against `now`.  The plate at `t_hit` is strictly between the two
    // frames and neither is it; `now` is the one the arrival was detected
    // against, and using the other would let a ball be found touching a plate
    // that had already moved on.
    const Eigen::Matrix<double, 6, 1> q = plateFrame(hit, now, ball_radius);

    // `q(5)` is the RELATIVE normal velocity at the contact point — `plateFrame`
    // subtracts the plate's own motion there before rotating into the plate's
    // axes — which is the quantity restitution acts on.  Negative is
    // approaching.
    const double approach = q(5);
    const double rebound = -kRestitution * approach;
    if (report) report->impact_approach = approach;

    // A bounce is a claim about the plate's velocity, exactly as a separation
    // is, so it is refused on the same grounds and for the same reason: near a
    // singularity `c_dot` and `omega` are arithmetic rather than physics, and a
    // rebound built on them would be a launch the mechanism never performed.
    // Landing is the conservative answer — it is what the ball did before this
    // model could bounce at all, and it cannot throw anything.
    if (now.rates_trustworthy && approach < 0.0 &&
        rebound >= bounceFloorSpeed(gravity, dt)) {
        // Reflect the normal component of the RELATIVE velocity and leave the
        // tangential part alone.  Adding `-(1 + e)` times the approach along
        // the normal does exactly that, and does it to the world velocity
        // without needing to decompose it: the plate's own contribution
        // cancels out of the difference.
        hit.flight_v -= (1.0 + kRestitution) * approach * now.normal();

        // Fly out the rest of the frame under the rebound.
        const double rest = dt - t_hit;
        hit.flight_p += hit.flight_v * rest + 0.5 * g * rest * rest;
        hit.flight_v += g * rest;
        return hit;
    }

    // The train is over: either the rebound is too small for this frame rate to
    // represent, or the ball arrived with no approach speed to reflect.  It
    // rolls on from where it touched down, which is the impact point rather
    // than wherever the rest of the frame would have carried it.
    out.airborne = false;
    out.rolling << q(0), q(1), q(3), q(4);   // z and vz are given up here
    return out;
}

}  // namespace caliburn
