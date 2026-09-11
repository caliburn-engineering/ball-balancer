// src/ball_contact.h
#pragma once

#include "ball_sim.h"
#include "table_kinematics.h"

#include <Eigen/Core>
#include <array>

namespace caliburn {

/// The plate's rigid-body motion in the world frame at one instant.
///
/// `RollingBallDynamics` needs none of this: it models a ball glued to a
/// surface, and a glued ball does not care how the surface moves.  A ball that
/// can come OFF needs all of it, because whether it stays is a question about
/// the surface's acceleration, not its tilt.
///
/// Everything here is analytic.  Nothing is differenced — the servo lag has a
/// closed-form derivative and the velocity Jacobian carries it through to the
/// pose, so the only finite difference in the whole contact calculation is the
/// single one in `normalAccel`, and it differences THESE rather than positions.
struct PlateMotion {
    Eigen::Matrix3d R = Eigen::Matrix3d::Identity();  ///< plate orientation
    Eigen::Vector3d c = Eigen::Vector3d::Zero();      ///< table centre, world
    Eigen::Vector3d c_dot = Eigen::Vector3d::Zero();
    Eigen::Vector3d omega = Eigen::Vector3d::Zero();  ///< angular velocity

    /// The accelerations, and they are the whole question — a ball leaves a
    /// plate because of how the surface ACCELERATES, not how it moves.
    ///
    /// **These are analytic, and that is not a refinement.**  They were a
    /// finite difference of `c_dot` and `omega` across consecutive frames,
    /// which is a defensible estimator for a smooth signal and these are not
    /// smooth: the leg rate is `(cmd - alpha) / tau`, and `cmd` is a zero-order
    /// hold that STEPS every time the loop changes its mind.  Differencing a
    /// step gives `1/dt`, so the estimator answered with a delta function
    /// wherever the command jumped.
    ///
    /// Measured at the corner of a square path on a thirty-second lap — a
    /// setpoint crawling at 24 mm/s, about as gentle as this demo gets — the
    /// difference reported `|omega_dot| = 101 rad/s^2` and a normal force of
    /// `-7.3 m/s^2`, and the ball was thrown off a plate that was barely
    /// moving.  A corner is a step in the reference velocity by construction
    /// (`pathVelocity` says so), so this fired on every corner of every lap.
    ///
    /// Analytically the servo lag differentiates in closed form: `cmd` is held
    /// across the frame, so `alpha_ddot = -alpha_dot / tau` exactly, and the
    /// same Jacobian that carries `alpha_dot` to `pose_dot` carries this to
    /// `pose_ddot`.  The one term left to difference is the Jacobian's own
    /// change, and that is a function of the leg angles and the pose — both
    /// continuous, neither of them stepping.  The corner then reads
    /// `|omega_dot| = 50 rad/s^2`, which is what a servo with a 0.05 s lag
    /// actually does when its command jumps.  See #23.
    Eigen::Vector3d c_ddot = Eigen::Vector3d::Zero();
    Eigen::Vector3d omega_dot = Eigen::Vector3d::Zero();

    /// The two maps this frame, kept so the NEXT frame can difference them.
    ///
    /// `J_v` takes leg rates to pose rates and `A` takes the tilt rates to the
    /// angular velocity.  Both depend only on the leg angles and the pose, so
    /// differencing them is differencing something continuous — which is the
    /// whole point of keeping them rather than differencing the rates.
    Eigen::Matrix3d J_v = Eigen::Matrix3d::Zero();
    Eigen::Matrix<double, 3, 2> A = Eigen::Matrix<double, 3, 2>::Zero();

    /// Whether `c_dot` and `omega` are worth believing.
    ///
    /// They come through the velocity Jacobian, which is `-J_pose^-1 J_alpha`,
    /// and near a kinematic singularity that inverse amplifies without bound:
    /// measured, an over-aggressive gain drives the plate to a condition
    /// number of 2464 and the rates that come back would launch the ball a
    /// metre into the air off a 0.26 m/s nudge.  The hop is arithmetic, not
    /// physics.
    ///
    /// So the rates carry their own credibility, and the contact test declines
    /// to act on rates that do not have it.  The threshold is the one the
    /// application already shows the user — condition number 20, the line
    /// where its own Jacobian readout turns red and says "Poor".
    bool rates_trustworthy = true;

    /// The plate's outward normal in world coordinates.
    Eigen::Vector3d normal() const { return R.col(2); }
};

// `kRatesUntrustworthyAbove` — the condition number above which `plateMotion`
// stops vouching for its own rates — moved to `table_kinematics.h` with #29.
// It was never only the contact model's: `maxBallAccel` sizes the corner
// fillet with it (#31) and `legCommand` bounds the plate's own travel with it,
// and a threshold about the velocity Jacobian belongs beside the velocity
// Jacobian.

/// Assemble the plate's motion from where the legs are and how they are moving.
///
/// `alpha_dot` and `alpha_ddot` are both exact rather than differenced: the
/// servos are a first-order lag, so `alpha_dot = (cmd - alpha) / tau` is the
/// model's own derivative and `alpha_ddot = -alpha_dot / tau` its own second,
/// for as long as `cmd` is held — which is the whole frame.  `servoAccel`
/// writes that down so no caller has to.
///
/// `prev` supplies the previous frame's `J_v` and `A` so their rates of change
/// can be differenced.  Passing `nullptr` drops those two terms, which is right
/// for the first frame — nothing has changed yet — and is why a plate assembled
/// without a history reports the acceleration its legs are producing rather
/// than nothing at all.
PlateMotion plateMotion(const TableKinematics& tk,
                        const TablePose& pose,
                        const std::array<double, 3>& alpha_rad,
                        const std::array<double, 3>& alpha_dot_rad_s,
                        const std::array<double, 3>& alpha_ddot_rad_s2 = {0.0, 0.0, 0.0},
                        const PlateMotion* prev = nullptr,
                        double dt = 0.0);

/// The servo lag's own second derivative: `alpha_ddot = -alpha_dot / tau`.
///
/// Exact while the command is held, which it is for the whole frame.  Here
/// rather than at each call site because five harnesses and the application all
/// need it and a sixth opinion about the servo model is how they come to
/// disagree.
std::array<double, 3> servoAccel(const std::array<double, 3>& alpha_dot_rad_s,
                                 double tau);

/// The normal force per unit mass holding the ball on the plate, `N / m`.
///
/// A ball is held on a surface by whatever normal force is needed to make its
/// centre follow the surface.  Write the centre's world position as
/// `P = c + R s`, with `s = (x, y, r)` the centre in the plate's own frame,
/// and take the component of Newton's second law along the plate normal:
///
///     N/m = g (n . z) + n . c_ddot
///           + n . [omega_dot x (R s) + omega x (omega x (R s))]
///           + 2 n . (omega x R s_dot)
///
/// The first two terms are the ones intuition supplies — gravity's share along
/// the normal, and the plate being driven up or down under the ball.  The
/// third is the plate's rotation swinging the contact point vertically, and
/// the fourth is Coriolis: a ball rolling on a tilting plate is held by a
/// different force than a ball sitting still on it.
///
/// **`N/m <= 0` is separation.**  The surface can push a ball, never pull it,
/// so a negative answer means the plate is accelerating away faster than
/// gravity can carry the ball after it, and contact is over.
///
/// `omega_dot` and `c_ddot` are read off `plate` rather than differenced out of
/// two frames.  See `PlateMotion::c_ddot` for why that stopped being a detail:
/// the leg command is a zero-order hold, so differencing the rates it drives
/// answers with a delta function at every step of the command.
double normalAccel(const PlateMotion& plate,
                   const Eigen::Vector3d& s,
                   const Eigen::Vector2d& s_dot,
                   double gravity);

/// The normal force per unit mass a plate would exert if it were holding still:
/// gravity's share along the normal, `g (n . z)`, and nothing else.
///
/// This is `normalAccel` with every rate term dropped, and it exists for the
/// one case where the rates cannot be believed.  A plate near a kinematic
/// singularity reports velocities that are arithmetic rather than physics, and
/// `stepBallContact` declines to act on them — but the ball is still resting on
/// something, and its rolling resistance still has to be scaled by something.
/// The pose is trustworthy even when the rates are not, so the tilt is the part
/// of the answer that survives.
double quasiStaticNormalAccel(const PlateMotion& plate, double gravity);

/// Where the ball is, and whether the plate is still touching it.
///
/// Two phases, two frames, and the frame is not an implementation detail.
/// While the ball rolls, its natural coordinates are the plate's — that is the
/// frame `RollingBallDynamics` integrates in and the frame the cascade plant's
/// state vector is written in, and neither should have to change because the
/// ball can now leave.  While it is airborne the only force on it is gravity,
/// which is a statement about the WORLD frame; re-expressing that in a frame
/// which is itself tilting and heaving would be a change of variables with
/// nothing whatever to recommend it.
///
/// So each phase keeps its own truth and `plateFrame` converts on demand, for
/// the things that want one answer regardless — the plots, the readout, and
/// the test of whether the ball is still over the plate.
struct BallState {
    bool airborne = false;

    /// [x, y, vx, vy] in the plate frame.  The truth while rolling.
    Eigen::Vector4d rolling = Eigen::Vector4d::Zero();

    /// Ball centre and velocity in the world.  The truth while airborne.
    Eigen::Vector3d flight_p = Eigen::Vector3d::Zero();
    Eigen::Vector3d flight_v = Eigen::Vector3d::Zero();
};

/// The ball's centre and velocity in the plate frame, in either phase:
/// `[x, y, z, vx, vy, vz]`.  While rolling, `z` is the ball's radius and `vz`
/// is zero — that is what rolling means.
Eigen::Matrix<double, 6, 1> plateFrame(const BallState& b,
                                       const PlateMotion& plate,
                                       double ball_radius);

/// The ball's centre and velocity in the world, in either phase.  This is what
/// a launch is made of: a ball leaving a moving plate carries the plate's
/// velocity at the contact point, not just its own rolling velocity.
void worldOf(const BallState& b, const PlateMotion& plate, double ball_radius,
             Eigen::Vector3d* p, Eigen::Vector3d* v);

/// Advance the ball one frame, changing phase when contact is lost or made.
///
/// Rolling while the plate can still push (`N/m > 0`), ballistic when it
/// cannot — and rolling regardless while `now.rates_trustworthy` is false,
/// because a separation is a claim about the plate's velocity and this is the
/// case where nobody knows what that is.  The launch takes the ball's full
/// world velocity — including the plate's motion at the contact point, which is
/// most of it when the legs are slamming.
///
/// The landing is inelastic: the normal component of the approach is absorbed
/// and the tangential part carries on rolling.
///
/// Inelastic rather than bouncing, deliberately.  A real ball does bounce, but
/// a restitution coefficient is a number nobody here has measured, and the
/// visible behaviour this ticket is about — the ball leaving the plate at all —
/// does not depend on it.  A bounce would be a second guess stacked on the
/// first.
///
/// `normal_accel_used`, if given, receives the `N/m` this frame was decided by
/// — the number, not a caller's reconstruction of it.  Which value that is
/// depends on the branch taken, and that is the point of reporting it:
///
///   - **In contact, rates believed:** `normalAccel` against the ball's state
///     at the frame's START, which is where the ball was when the plate was
///     either still pressing it or not.  Positive held contact; negative ended
///     it, and is the frame the ball left on.
///   - **In contact, rates not believed:** `quasiStaticNormalAccel`, which is
///     what actually scaled the rolling resistance — the tilt is trustworthy
///     when the rates are not, so it is the part of the answer that survives.
///   - **Already in flight:** zero.  A plate touching nothing presses with
///     nothing.
///
/// Anything measuring the margin the plate is holding the ball by has to ask
/// for it here.  Recomputing `normalAccel` outside is a second opinion about
/// which state and which guard it was evaluated under, and it will differ:
/// against the end-of-frame ball rather than the start-of-frame one, without
/// the trust gate, and with a `z` that has to be assumed.
BallState stepBallContact(const RollingBallDynamics& dynamics,
                          const BallState& b,
                          const PlateMotion& now,
                          const PlateMotion& prev,
                          const TablePose& pose,
                          double ball_radius,
                          double gravity,
                          double dt,
                          double* normal_accel_used = nullptr);

}  // namespace caliburn
