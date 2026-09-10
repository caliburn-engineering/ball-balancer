// src/auto_balance.h
#pragma once

#include "table_kinematics.h"

#include <Eigen/Core>
#include <array>
#include <cmath>

namespace caliburn {

/// A designed gain, plus everything about the operating point it was designed
/// about.  The plate simulates a *nonlinear* mechanism; the gain came from a
/// linearisation.  Every field here is part of the contract between them, and
/// the loop is only honest while all of them agree.
struct AutoBalanceDesign {
    /// 3 x 7 state feedback, u = -K(x - x_ref), u in leg-command deviations.
    Eigen::MatrixXd K;

    /// The leg angle the plant was linearised about — the `a0` parameter of
    /// the cascade preset.  States 0..2 are deviations *from this*, so a
    /// mismatch does not fail, it quietly regulates to the wrong pose.
    double home_leg_rad = M_PI / 4.0;

    /// First-order leg lag, the `tau` of the same preset.  The simulator had
    /// no servo dynamics at all before the loop was closed: the slider WAS the
    /// leg angle.  Closing u = -Kx around that is a unit-delay algebraic loop
    /// on the leg states, so the lag the model claims has to actually exist.
    double servo_tau = 0.05;

    /// Servo travel.  A command outside it is clamped, not refused — the plate
    /// really does stop there, and hiding that would make an over-aggressive
    /// tuning look better than it is.
    double alpha_min_rad = 10.0 * M_PI / 180.0;
    double alpha_max_rad = 80.0 * M_PI / 180.0;

    /// The plant the gain was designed against.  The model panel's physical
    /// sliders move it and the simulated plate does not follow, so a gain
    /// designed for longer legs — or for the moon — is not wrong-shaped, it is
    /// simply wrong, and nothing downstream could tell.  Default-constructed
    /// to zeros, so a design nobody filled in is refused rather than trusted.
    TableParams mechanism{};
    double gravity = 0.0;  ///< [m/s^2]; zero is "unset", never a real plant
};

/// Do two designs describe the same plate?  Geometry and gravity: the shape
/// the legs make and the acceleration acting on the ball are independent, and
/// getting either wrong is the same silent failure.  `alpha_min`/`alpha_max`
/// are travel limits rather than shape, and are not compared.
bool samePlant(const TableParams& a, double g_a,
               const TableParams& b, double g_b);

/// True when `d.K` has the shape the cascade plant's state feedback must have.
/// Shape is necessary, not sufficient — the caller also has to know the plant
/// IS the cascade, which no matrix dimension can tell it.
bool gainFitsCascade(const AutoBalanceDesign& d);

/// Assemble the cascade plant's state vector out of what the simulator holds:
///
///     [da1, da2, da3, x, y, x', y']
///
/// Leg deviations from the linearisation point first, then the ball's position
/// and velocity in the plate's own frame — the frame `RollingBallDynamics`
/// already integrates in, so nothing is transformed here.
Eigen::VectorXd cascadeState(const std::array<double, 3>& alpha_rad,
                             double home_leg_rad,
                             const Eigen::Vector4d& ball);

struct LegCommand {
    std::array<double, 3> alpha_rad;
    bool saturated;  ///< at least one leg hit its travel limit

    /// The command asked for a pose the mechanism cannot assemble, and was
    /// scaled back until it could.  Distinct from `saturated`: a triple can be
    /// inside every servo's travel and still not exist, because the travel
    /// limits are a box and the workspace is not.
    bool clipped_to_workspace;
};

/// What the loop is asked to make the ball do: be somewhere, and be moving.
///
/// The position half is the setpoint the loop has always had.  The velocity
/// half is the **reference velocity**, and it is zero for a setpoint being
/// held — which is why the loop went so long without one.
///
/// A ball at rest anywhere on a flat plate is an equilibrium of this plant, so
/// a STATIONARY reference state is reachable with zero leg deviation and the
/// regulator tracks it with no feedforward at all.  That stops being true the
/// moment the setpoint moves: `[.., x_sp, y_sp, 0, 0]` then claims the ball
/// should be at the setpoint *and stationary*, which is false, and the loop
/// spends its effort fighting the very motion it was asked for.  Measured on a
/// 120 mm circle at a ten-second lap: **3.66 mm of mean error with the
/// reference velocity supplied, 35.52 mm without** — nearly ten times.  See
/// #24 and `setpoint_path.h`.
///
/// **It has no gain of its own, and must not be given one.**  It enters as a
/// reference velocity and is multiplied by K's velocity columns, which LQR
/// designs from the `x' ball` and `y' ball` weights — so it is already tunable,
/// from the two sliders that own those numbers.  A separate feedforward gain
/// would be a second opinion about a number K owns.
struct BallReference {
    Eigen::Vector2d position{Eigen::Vector2d::Zero()};   ///< [m], plate frame
    Eigen::Vector2d velocity{Eigen::Vector2d::Zero()};   ///< [m/s], plate frame
};

/// The loop, evaluated once: u = -K(x - x_ref), commanded as `home + u`.
///
/// `x_ref` is `[0, 0, 0, ref.position, ref.velocity]` — the whole reference
/// state, assembled here rather than by the caller.  It lived in the caller
/// once, as a bias on the *measurement* handed in, and the arithmetic was
/// identical because only the difference enters `u = -K(x - x_ref)`.  The
/// trouble was that the loop's own header then said it needed no feedforward
/// while the application supplied one two files away, and the test harness
/// carried a third copy.  See `BallReference` and #24.
///
/// A gain of the wrong shape returns the home pose and no saturation, so a
/// caller that forgets to check `gainFitsCascade` gets a flat plate rather
/// than an out-of-bounds read.
///
/// `tk` is here because clamping each leg to its own travel is not enough.
/// Barely half of the servo box has an assembly at all, so a per-leg clamp can
/// hand back a triple the plate cannot make — and a plate that cannot make its
/// command has no pose, so it freezes at whatever tilt it last held and rolls
/// the ball off.  That was the second half of issue #22.
///
/// The command is scaled back toward the home pose until it is assemblable,
/// which gives up magnitude and keeps direction — what saturation ought to do.
/// The level pose always assembles, so the search always has an answer.
LegCommand legCommand(const TableKinematics& tk,
                      const AutoBalanceDesign& d,
                      const std::array<double, 3>& alpha_rad,
                      const Eigen::Vector4d& ball,
                      const BallReference& ref);

/// Pull `target` back toward `safe` until the mechanism can be assembled there.
///
/// The servo travel limits are a box and the workspace is not — barely half
/// the box has an assembly at all, and it is under no obligation to be convex.
/// So both a command and a servo step can land outside it, and a leg triple
/// with no assembly gives the plate no pose: it freezes at whatever tilt it
/// last held, and a frozen tilted plate rolls the ball off.  That was issue
/// #22's second half.
///
/// `safe` must itself be assemblable, and every caller has one to hand: the
/// home pose for a command, the legs' present position for a step.  The legs
/// start at home and are never moved anywhere unassemblable, so that second
/// one holds by induction.
///
/// Scaling gives up magnitude and keeps direction, which is what saturation
/// ought to do.  Bisection, but the answer does not depend on the workspace
/// being convex along the ray: the returned point is one that was actually
/// tested, never an interpolated boundary.
struct Retreat {
    std::array<double, 3> alpha_rad;
    bool retreated;  ///< the target had no assembly and was pulled back

    /// `safe` did not assemble either, so there was nothing to retreat TO and
    /// `alpha_rad` is just `safe` handed back.  The precondition is broken and
    /// the plate is about to freeze; this is the flag that says so rather than
    /// letting it look like an ordinary clip.  It cannot happen through the
    /// two call sites here — the level pose always assembles, and the legs are
    /// never moved anywhere that does not — but "cannot happen" is worth one
    /// bool when the alternative is a silent stall.
    bool safe_was_unassemblable;
};

Retreat retreatToWorkspace(const TableKinematics& tk,
                           const std::array<double, 3>& target,
                           const std::array<double, 3>& safe);

/// Where the ball will be when it comes back down, in the plate's frame.
///
/// A ball in the air is one the plate cannot touch, so `u = -Kx` on where it
/// IS regulates a quantity nothing can move.  Where it will LAND is a
/// different matter: the plate has the whole flight to get underneath it, and
/// the vertical state is exactly what says how long that is.
///
/// Solves `z + vz t - g t^2 / 2 = r` for the positive root and carries the
/// horizontal motion forward.  Returns the present position unchanged for a
/// ball already at or below the surface, because a prediction nobody can act on
/// is worse than none.
///
/// A ball above the surface always lands: with `dz > 0` the discriminant
/// `vz^2 + 2 g dz` is positive however hard the ball was thrown, so there is no
/// such thing here as escaping the plate upward.  The implementation still
/// tests it, as an assertion rather than a case.
///
/// Approximate on purpose: it ignores the plate's own motion during the
/// flight, which is real but which the plate is about to choose.  See #23.
Eigen::Vector2d predictedLanding(const Eigen::Matrix<double, 6, 1>& ball_plate,
                                 double ball_radius,
                                 double gravity);

/// Advance the three first-order servos one step, integrated exactly:
///
///     alpha <- cmd + (alpha - cmd) * exp(-dt / tau)
///
/// Exact rather than forward-Euler because the plate runs at a fixed 60 Hz
/// against tau = 0.05 s — three steps per time constant, where Euler already
/// overshoots, and where a user dragging tau below 1/30 s would make Euler
/// oscillate and then diverge.  The exponential form cannot, at any dt.
std::array<double, 3> stepServos(const std::array<double, 3>& alpha_rad,
                                 const std::array<double, 3>& cmd_rad,
                                 double tau,
                                 double dt);

/// The same step, stopped at the workspace boundary.
///
/// `stepServos` lags each leg independently, so its result sits on the
/// straight line from where the legs are to where they were told to go — and
/// that line can leave the workspace even when both of its ends are inside.
/// Clipping the COMMAND is therefore not enough, which is what the last two
/// failing kick directions turned out to be: one frame, mid-flight, with no
/// assembly.
///
/// The legs stop where the mechanism stops them, which is what a real one
/// does when it is driven into a bind.
std::array<double, 3> stepServosOnPlate(const TableKinematics& tk,
                                        const std::array<double, 3>& alpha_rad,
                                        const std::array<double, 3>& cmd_rad,
                                        double tau,
                                        double dt);

/// The weights the LQR designer opens on.
///
/// Unit weights are the textbook default, and on this plant they do not
/// balance the ball.  The rolling model carries Coulomb resistance, so a plate
/// tilted by less than atan(c_rr) = 0.573 deg cannot start the ball moving at
/// all — and a state feedback with no integral term has no way out of that
/// dead band.  It parks the ball wherever the tilt it is asking for falls
/// inside it, which under unit weights is 37 mm off centre: a loop that looks
/// broken while behaving exactly as designed.
///
/// The residual is inversely proportional to the position gain, so the fix is
/// a default that actually weights ball position — 1 mm under these.  Pinned
/// by `test_auto_balance`, both the good case and the dead band itself.
///
/// Keyed on the plant's dimensions, the only thing this layer can see.  A
/// different 7-state, 3-input plant gets the ball-balancer's opening tuning,
/// which is a poor guess rather than a wrong answer — the sliders are right
/// there, and every other shape still opens on unit weights.
Eigen::VectorXd defaultLqrStateWeights(int n);
Eigen::VectorXd defaultLqrInputWeights(int m);

/// The largest disturbance the interface can hand the loop, in m/s, as a SPEED
/// in any direction.
///
/// Here rather than beside the slider because it is what bounds an honest
/// preset.  "This tuning does not lose the ball" is only a claim if something
/// says how hard the ball can be hit, and the answer is whatever the UI lets a
/// visitor ask for.  The slider reads this constant and so does the test that
/// checks the presets against it — raise it and the test fails, which is the
/// point.
///
/// **It was 0.5, and 0.5 was measured against the wrong situation.**  That
/// number came from shoving a ball sitting at rest at the centre of a level
/// plate with a held setpoint, and from rest every shipped tuning does survive
/// it.  The demo does not open like that: it opens tracking a circle, so the
/// legs are already displaced and the ball is already running at 75 mm/s when
/// the visitor reaches for Nudge.  A shove there lands on a plate that is
/// already working, and the envelope is much smaller — measured over 36 points
/// of the lap and 24 directions, every preset is clean at 0.30 and Nominal
/// loses the ball at 0.35.  See
/// `test_a_shove_while_tracking_is_rejected_from_every_direction`.
///
/// The slivers thin as the shove shrinks rather than stopping at a threshold —
/// the same shape as `kMaxSetpointSpeed`'s bound, and the same reason for
/// leaving margin rather than sitting on the first clean measurement.
constexpr double kMaxNudgeSpeed = 0.30;

/// What one Nudge button may add along its own axis.
///
/// **Derived, not chosen, because the two buttons compose.**  "Nudge +x" and
/// "Nudge +y" each add their slider's worth to one axis, so pressing both gives
/// a shove of `sqrt(2)` times the slider at 45 degrees — and the slider's top
/// used to BE `kMaxNudgeSpeed`, so two clicks handed the loop 0.707 m/s against
/// a constant whose own comment called itself the largest disturbance the
/// interface can hand it.  Measured, that pair lost the ball at every one of 36
/// points of the lap.
///
/// Dividing here is what makes that sentence true: the worst the pair can
/// compose to is exactly `kMaxNudgeSpeed`, which is the number the envelope was
/// measured at.  Bounding the buttons rather than clamping the ball afterwards,
/// because a clamp would also be silently deciding what the ball's own tracking
/// velocity is allowed to be, and that is not a disturbance.
inline const double kMaxNudgePerAxis = kMaxNudgeSpeed / M_SQRT2;

/// A named tuning, as a visitor meets it.
///
/// The demo's whole argument is that controller design has consequences you
/// can watch, and the fastest way to make that argument is to let someone flip
/// between two tunings of the same plant.  A weighting matrix cannot make it
/// in five seconds; a pair of buttons can.
///
/// A *tuning*, not a plant: "preset" means a built-in plant model everywhere
/// else in this application, and these choose Q and R against one plant rather
/// than choosing the plant.
struct LqrPreset {
    const char* name;

    /// One line, in the UI, saying what is about to happen.  Written from the
    /// measurements in `test_auto_balance`, so a claim here that stops being
    /// true fails a test rather than merely misleading somebody.
    const char* blurb;

    double q_position;  ///< Q on ball x and ball y
    double q_velocity;  ///< Q on x' and y'
    double r;           ///< R on every leg command
};

/// The three tunings, in the order the UI offers them: sluggish, shipped,
/// saturating.  Ordered rather than sorted by any number — it is the order a
/// visitor should press them in.
///
/// Every value was measured against the nonlinear plate rather than chosen for
/// roundness; the measurements and what they rule out are recorded beside the
/// table in `auto_balance.cpp`, and the claims each tuning makes are asserted
/// in `test_auto_balance`.
const std::array<LqrPreset, 3>& lqrPresets();

/// The tuning the application opens on.  `defaultLqrStateWeights` is defined
/// in terms of this rather than the other way round, so "Nominal is the
/// startup tuning" is arithmetic rather than a coincidence two literals have
/// to keep agreeing about.
const LqrPreset& nominalPreset();

/// A preset's weights, sized to a plant.
///
/// Keyed on dimensions, exactly as `defaultLqrStateWeights` is and for the
/// same reason: this layer can see the plant's shape and nothing else.
/// Anything that is not 7 states gets unit weights — a poor guess rather than
/// a wrong answer, and the UI does not offer the presets off the cascade
/// anyway, which it establishes by NAME (see `isCascadeModel`).
Eigen::VectorXd presetStateWeights(const LqrPreset& p, int n);
Eigen::VectorXd presetInputWeights(const LqrPreset& p, int m);

/// Which preset these weights ARE, or -1 for none.
///
/// The presets are an offer, not a mode: the weights stay editable after one
/// is chosen, and the moment a slider moves the answer is -1 again.  There is
/// therefore no "selected preset" to store — storing one would go stale on the
/// first drag and the UI would keep claiming a tuning the plant no longer has.
/// This asks the weights instead, every frame, which cannot disagree with them.
///
/// Lives here rather than in the panel because it is the rule the highlight
/// means, and a rule worth stating is worth testing.
///
/// -1 for anything that is not 7 states by 3 inputs: at any other shape every
/// preset is unit weights, so a "match" would only be reporting that they have
/// stopped differing.  Shape is necessary and not sufficient, exactly as in
/// `gainFitsCascade` — the caller also has to know the plant IS the cascade,
/// which no vector length can tell it.
int activePreset(const Eigen::VectorXd& q, const Eigen::VectorXd& r);

}  // namespace caliburn
