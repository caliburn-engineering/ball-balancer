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
    /// smooth: the leg rate is `servoRate`, which away from its limit is
    /// `(cmd - alpha) / tau`, and `cmd` is a zero-order hold that STEPS every
    /// time the loop changes its mind.  Differencing a
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
    /// across the frame, so `alpha_ddot = -alpha_dot / tau` exactly while the
    /// lag is what is limiting the leg, and the same Jacobian that carries
    /// `alpha_dot` to `pose_dot` carries this to `pose_ddot`.  The one term
    /// left to difference is the Jacobian's own change, and that is a function
    /// of the leg angles and the pose — both continuous, neither of them
    /// stepping.  The corner then reads `|omega_dot| = 50 rad/s^2`, which is
    /// what a servo with a 0.05 s lag actually does when its command jumps.
    /// See #23.
    ///
    /// **The lag is not always what is limiting the leg.**  Since #32 the servo
    /// is rate-limited too, and while it is ramping `alpha_dot` is a constant
    /// so `alpha_ddot` is exactly zero — these two vectors go to zero with it,
    /// in the frames the loop is slamming hardest.  `servoRate` and
    /// `servoAccel` are the single statement of both derivatives; a caller
    /// deriving either from `(cmd - alpha) / tau` itself gets the unlimited
    /// servo back and does not notice.
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
/// servos are a rate-limited first-order lag, so `servoRate` is the model's own
/// derivative and `servoAccel` its own second, for as long as `cmd` is held —
/// which is the whole frame.  The two write that down so no caller has to.
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

/// The rate-limited servo's own first derivative:
///
///     alpha_dot = clamp((cmd - alpha) / tau, -rate_max, +rate_max)
///
/// Exact while the command is held, which it is for the whole frame.  Here
/// rather than at each call site because five harnesses and the application all
/// need it and a sixth opinion about the servo model is how they come to
/// disagree — and the clamp above is the same handover `stepServos` integrates
/// through, so the rate the plate reports is the rate the legs are actually
/// moving at.
///
/// A lagless servo (`tau <= 0`) is driven entirely by the limit: it closes any
/// error at `rate_max` and reports zero only once there is none.  Without a
/// limit its rate is unbounded and so unrepresentable, and it reports zero, as
/// it did before the limit existed.
std::array<double, 3> servoRate(const std::array<double, 3>& alpha_rad,
                                const std::array<double, 3>& cmd_rad,
                                double tau,
                                double rate_max);

/// The rate-limited servo's own second derivative: `-alpha_dot / tau` while the
/// lag is in charge, and **exactly zero while the rate limit is**.
///
/// That zero is the whole reason the rate limit is worth having.  `c_ddot` and
/// `omega_dot` are the terms that decide whether the ball separates, and they
/// are driven by `alpha_ddot` — which now vanishes precisely when the loop is
/// slamming the legs hardest.  At the handover `alpha_ddot` jumps to
/// `-alpha_dot / tau`, the same magnitude an unlimited servo would have had,
/// but later and from a smaller error: strictly better, never worse.  See #32.
///
/// Saturation is read off `alpha_dot` rather than passed in, which is exact
/// because `servoRate` clamps to `rate_max` itself and so lands on it bit for
/// bit.  That clamp is also why the test is `|alpha_dot| >= rate_max` rather
/// than `>`: a ramping leg reports exactly the limit, so `>` would find no
/// saturated leg at all.
///
/// At the handover exactly — `|cmd - alpha| = rate_max * tau`, where the lag
/// alone would produce precisely `rate_max` — the two readings coincide and
/// this reports the ramp's zero.  A measure-zero tie, broken toward the side
/// the leg is arriving FROM; `stepServos` puts that same instant on the decay
/// branch, so the pair disagree at one point and nowhere else.
std::array<double, 3> servoAccel(const std::array<double, 3>& alpha_dot_rad_s,
                                 double tau,
                                 double rate_max);

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

/// How fast the plate's contact point under the ball is RISING, in m/s along
/// the plate normal: `n . (c_dot + omega x R s)`.
///
/// **This is the `p` of `u_rebound = e u + p (1 + e)`, and its sign is the
/// whole of whether the loop pumps a bouncing ball or damps one.**  A plate
/// rising into a falling ball hands it back more than it arrived with, at
/// `1.94x` the plate's own speed at `e = 0.94`; a plate that is level or
/// descending at the moment of impact cannot, and the apex ratio is then at
/// most `e^2` whatever the loop is doing.  It is a hard result rather than a
/// tuning, which is why it is worth a named function and a report field
/// instead of being buried inside the bounce.
///
/// `s` is the ball's CENTRE in the plate frame, not the point of the surface
/// beneath it, and the two give the same answer: they differ by `r n`, and
/// `n . (omega x r n)` is zero.  Taking the centre is what makes this the same
/// arithmetic `plateFrame` already subtracts out of an airborne ball's
/// velocity, rather than a second opinion about it.
double contactNormalRate(const PlateMotion& plate, const Eigen::Vector3d& s);

/// The coefficient of restitution for this ball arriving on this plate.
///
/// **Cited, not chosen.**  Chai, Y. et al., "Restitution coefficient of various
/// particles based on acoustic technology", *J. Phys.: Conf. Ser.* **2557**
/// 012057 (2023) — a sphere dropped down a guide onto a 304 stainless-steel
/// base plate, with the intervals between successive impacts timed acoustically
/// at 48 kHz.  Their POM sphere gives **e = 0.94**.
///
/// It applies to *this* ball because the ball already had a material, implicitly
/// and unavoidably: `kPlateBall` fixes 50 g at a 20 mm radius, which is
/// `rho = 1492 kg/m^3` — the engineering-polymer range, and emphatically not
/// steel, aluminium or a hollow shell.  Chai's POM sphere is 1410 kg/m^3.  So
/// the density argument and the coefficient come from the same place, and
/// nothing about the ball had to be retuned to accept the number.
///
/// **Two extrapolations, recorded rather than glossed.**  Their drop is 200 mm,
/// so their first impact is about 1.98 m/s, where the separations here are
/// nearer 0.2 m/s; restitution is velocity-dependent and generally rises as
/// impact speed falls, so 0.94 is more likely a floor than a ceiling at these
/// speeds (Schwager & Poeschel, arXiv:1204.0001, record King et al. finding the
/// dependence explicitly non-monotonic at low speed for POM on steel).  And
/// their sphere is 14.8 mm across against this one's 40 mm — the same order,
/// not the same ball.
///
/// It lives beside the contact model rather than in `BallParams` because
/// restitution is a property of the *pair* of materials that meet, not of the
/// ball on its own.  A different plate would change it without the ball
/// changing at all.
inline constexpr double kRestitution = 0.94;

/// The rebound speed below which a bounce is not resolvable at this frame rate,
/// and the ball is declared landed: `u < g dt / 2`.
///
/// **The termination rule is the hard part, not the coefficient.**  Repeated
/// restitution gives infinitely many bounces in finite time, so something has
/// to end the train — and the obvious rules are all height or velocity
/// thresholds picked against the passive bouncing case, which is exactly what
/// would kill the small deliberate hops a hopping controller commands (#34).
///
/// So the rule is not a chosen threshold at all.  A ball rebounding at `u` is
/// airborne for `2u/g`; below `u = g dt / 2` that whole flight fits inside one
/// frame, so the next sample already finds the ball back on the plate and the
/// integrator has nothing to integrate.  Such a bounce is not *small*, it is
/// *unrepresentable* — the simulation cannot show it whatever it decides.
///
/// Two things follow that a chosen threshold would not give.  It scales with
/// the frame rate rather than with the ball, so a faster simulation resolves
/// finer bounces instead of inheriting a number tuned at 60 Hz.  And at 60 Hz
/// it is 0.082 m/s — a 0.34 mm hop, an order below the millimetres the shipped
/// tuning's own passive hop reaches and further still below anything a
/// deliberate hop would ask for.  A ball that is genuinely hopping cannot reach
/// it.
inline constexpr double bounceFloorSpeed(double gravity, double dt) {
    return 0.5 * gravity * dt;
}

/// Advance the ball one frame, changing phase when contact is lost or made.
///
/// Rolling while the plate can still push (`N/m > 0`), ballistic when it
/// cannot — and rolling regardless while `now.rates_trustworthy` is false,
/// because a separation is a claim about the plate's velocity and this is the
/// case where nobody knows what that is.  The launch takes the ball's full
/// world velocity — including the plate's motion at the contact point, which is
/// most of it when the legs are slamming.
///
/// **Arrival is a bounce, at `kRestitution`.**  It is applied to the *relative*
/// normal velocity at the contact point rather than to the ball's world
/// velocity, which is the difference between a model and a decoration: a plate
/// rising into a falling ball throws it harder than it arrived, and a plate
/// running away softens the landing to nothing.  That same term is what lets a
/// controller hop the ball on purpose by driving the legs down and then up, and
/// it is what `holdContactDown` exists to keep the loop from doing by accident.
///
/// **This reverses a call #23 made on purpose, and the words being reversed
/// are these:** "Inelastic rather than bouncing, deliberately.  A real ball does
/// bounce, but a restitution coefficient is a number nobody here has measured,
/// and the visible behaviour this ticket is about — the ball leaving the plate
/// at all — does not depend on it."  The second clause still stands; the first
/// stopped being true when the coefficient was measured by somebody else and
/// cited.  A ball that arrives and sticks reads as a bug in the physics rather
/// than as a simplification of it, and that is a cost to a demo whose whole
/// argument is that you can watch the model be right.
///
/// **Tangential velocity is carried across the impact unchanged, deliberately.**
/// A tangential impulse is `mu` times the normal one at most, and its sign and
/// size depend on the ball's spin at the moment it lands — but this ball has no
/// spin state to consult.  `RollingBallDynamics` assumes rolling without
/// slipping, which fixes spin from velocity while the ball is *on* the plate
/// and says nothing whatever about a ball in flight.  Inventing a tangential
/// law would therefore be inventing the spin it depends on.  Frictionless in
/// the tangential direction is the one choice that assumes nothing — and for a
/// ball that left the plate rolling without slipping it is not even an
/// assumption, since it keeps that spin through the flight and arrives with a
/// contact point already at rest.
///
/// The train ends when a rebound falls below `bounceFloorSpeed`, at which point
/// the normal motion is given up and the ball rolls on from where it touched
/// down.  See that constant for why the floor is derived rather than chosen.
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
///
/// `impact_approach`, if given, receives the RELATIVE normal velocity the ball
/// arrived at if this frame resolved an impact, and zero if it did not.
/// Negative is approaching, so a bounce always reports a negative number.
///
/// **This is the quantity the no-pumping property is stated in.**
/// `u_rebound = e u + p (1 + e)`, so `p <= 0` gives `u_rebound <= e u` and
/// successive arrivals within one flight fall by at least `e` — which is the
/// apex ratio `e^2` in the form that survives a plate that is itself moving.
/// A plate-frame apex is not that quantity: a plate descending under a ball
/// makes the gap grow with no energy going into the ball at all.  Reported
/// rather than reconstructed for `normal_accel_used`'s reason — it is read at
/// the sub-frame crossing, which is an instant no caller has.
BallState stepBallContact(const RollingBallDynamics& dynamics,
                          const BallState& b,
                          const PlateMotion& now,
                          const PlateMotion& prev,
                          const TablePose& pose,
                          double ball_radius,
                          double gravity,
                          double dt,
                          double* normal_accel_used = nullptr,
                          double* impact_approach = nullptr);

}  // namespace caliburn
