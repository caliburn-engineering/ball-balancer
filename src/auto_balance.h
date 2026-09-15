// src/auto_balance.h
#pragma once

#include "ball_contact.h"
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

    /// The command asked for a pose the plate cannot be held at, and was
    /// scaled back until it could.  Distinct from `saturated`: a triple can be
    /// inside every servo's travel and still not exist, because the travel
    /// limits are a box and the workspace is not — and a triple can exist and
    /// still not be somewhere to steer to, because a second assembly waits on
    /// the far side of a singularity (#29).
    ///
    /// **`clipped_to_holdable`, not `clipped_to_workspace`**, which is what
    /// this was called until #29, for the same reason `retreatToHoldable` is
    /// not `retreatToWorkspace`: the workspace is the set that assembles, and
    /// the clip is against a smaller one.
    bool clipped_to_holdable;
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
/// The command is scaled back toward the home pose until the plate can be HELD
/// there, which gives up magnitude and keeps direction — what saturation ought
/// to do.  Held rather than merely assembled: a triple just past a singularity
/// has an assembly and it is the WRONG one, mirrored about the crossing, and a
/// plate steered onto it tilts away from the ball while the loop asks for the
/// opposite.  See `retreatToHoldable` and #29.  The level pose always
/// assembles and sits at a condition number of 4.7, so the search always has an
/// answer.
LegCommand legCommand(const TableKinematics& tk,
                      const AutoBalanceDesign& d,
                      const std::array<double, 3>& alpha_rad,
                      const Eigen::Vector4d& ball,
                      const BallReference& ref);

/// Pull `target` back toward `safe` until the plate can be held there.
///
/// **`retreatToHoldable`, not `retreatToWorkspace`**, which is what this was
/// called until #29.  The workspace is the set of triples that assemble, and
/// retreating into it is not enough — the last three per cent of it, against a
/// singularity, is where the plate comes back mirrored.  The old name would
/// now be describing the wrong set.
///
/// The servo travel limits are a box and the workspace is not — barely half
/// the box has an assembly at all, and it is under no obligation to be convex.
/// So both a command and a servo step can land outside it, and a leg triple
/// with no assembly gives the plate no pose: it freezes at whatever tilt it
/// last held, and a frozen tilted plate rolls the ball off.  That was issue
/// #22's second half.
///
/// **Existence is not the whole bound, and #29 is what the rest of it costs.**
/// A triple can have a perfectly good assembly and still be somewhere no plate
/// should be steered: approaching a direct-kinematics singularity, the built
/// assembly the plate is in meets a second one and they swap, so a pose
/// marched through the meeting comes out mirrored — tilting away from the ball
/// while the loop asks for the opposite.  Measured under the aggressive tuning,
/// the plate crossed a Jacobian condition number of 290 and spent the rest of
/// the run tilted downhill toward a ball it was trying to catch.  So the set
/// retreated into is `can_hold`'s rather than `can_assemble`'s, bounded by
/// `condition_limit`.
///
/// `safe` must itself be holdable, and every caller has one to hand: the home
/// pose for a command, the legs' present position for a step.  The legs start
/// at home — condition number 4.7 on the shipped plate — and are never moved
/// anywhere unholdable, so that second one holds by induction.
///
/// Scaling gives up magnitude and keeps direction, which is what saturation
/// ought to do.
///
/// **March, then bisect inside the bracket**, the same shape
/// `max_conditioned_tilt` uses.  Halving the whole ray would assume the
/// holdable points are one unbroken run from `safe`, and they are not — the
/// condition number rises toward a singularity and falls again past it, so a
/// ray can run good, bad, good, and a bisection that lands in the third band
/// returns a command the legs cannot travel to.  What comes back is a point
/// that was actually tested AND that the march reached by an unbroken run,
/// never an interpolated boundary.
struct Retreat {
    std::array<double, 3> alpha_rad;
    bool retreated;  ///< the target could not be held and was pulled back

    /// `safe` could not be held either, so there was nothing to retreat TO and
    /// `alpha_rad` is just `safe` handed back.  The precondition is broken and
    /// the plate is about to freeze; this is the flag that says so rather than
    /// letting it look like an ordinary clip.  It cannot happen through the
    /// two call sites here — the level pose always assembles and is far from
    /// any singularity, and the legs are never moved anywhere that is not —
    /// but "cannot happen" is worth one bool when the alternative is a silent
    /// stall.
    bool safe_was_unholdable;
};

/// `retreatToHoldable` with the safe end every leg COMMAND has: the level pose.
///
/// Two call sites want exactly this — the loop's own command and the
/// constrained one `holdContactDown` hands back — and they wanted it in the
/// same four lines, which is four lines of "which safe end, and which condition
/// limit" stated twice.  The level pose always assembles and sits at a
/// condition number of 4.7, so the search always has an answer; the limit is
/// `kRatesUntrustworthyAbove`, the same line the application's own Jacobian
/// readout turns red at.
std::array<double, 3> commandRetreatedToHoldable(const TableKinematics& tk,
                                                 const AutoBalanceDesign& d,
                                                 const std::array<double, 3>& cmd_rad,
                                                 bool* retreated = nullptr);

Retreat retreatToHoldable(const TableKinematics& tk,
                          const std::array<double, 3>& target,
                          const std::array<double, 3>& safe,
                          double condition_limit);

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

/// The leg command with its heave held down until the plate's contact point
/// under the ball cannot rise.  **The no-pumping constraint, and it is a proof
/// rather than a measurement.**
///
/// Restitution acts on the relative normal velocity at the contact point, so a
/// ball arriving at `u` onto a contact point rising at `p` leaves at
///
///     u_rebound = e u + p (1 + e)
///
/// `p <= 0` therefore gives `u_rebound <= e u` and an apex ratio of at most
/// `e^2`, **always**, whatever else the loop is doing.  A rising plate always
/// adds energy, at `1.94x` its own speed at `e = 0.94`.  That is a hard result
/// about the impact law and not a tuning, which is why the loop is given the
/// one scalar to constrain rather than a fourth target to chase.
///
/// **Why a constraint and not a better setpoint.**  Four airborne control laws
/// were measured against the 72-direction and 24-setting sweeps and all four
/// were worse than or equal to doing nothing: `predictedLanding` with the
/// ball's velocity zeroed (what ships), the landing point plus the real
/// velocity, `predictedRest`, and keeping the path feedforward through the
/// flight.  They were not four control laws.  `K`'s output is a leg TRIPLE,
/// which is simultaneously tilt (differential) and heave (common-mode), and
/// those do different jobs while the ball is in the air — tilt aims the normal
/// impulse at landing, heave sets its magnitude.  Feeding a single target
/// through `K` leaves their relative phase uncontrolled, so all four were
/// different ways of saying nothing about `p`.  Worst of them, `predictedLanding`
/// alone: the prediction collapses onto the ball at every impact and springs
/// out again at every rebound, so a bouncing ball hands the loop a target
/// oscillating at the bounce frequency, the plate rises about 0.02 m/s into
/// each arrival, and that feeds back roughly 5% per impact against the 12%
/// `e^2` takes out.  See the decision record's D8.
///
/// **Tilt authority is untouched, and exactly so.**  The correction is the leg
/// rate that `J_v` maps to pure heave — `J_v^-1 (0, 0, dz)` — so `phi_dot` and
/// `theta_dot` come out bit for bit unchanged and `omega` with them.  A uniform
/// `(1, 1, 1)` would have been the obvious reading of "common-mode" and is not
/// the same thing away from the symmetric pose, where it tilts the plate
/// slightly as it lifts it.
///
/// `p` is HOMOGENEOUS in the leg rates, which is what makes it enough to
/// evaluate this at the frame's open: the servo carries `cmd - alpha` across
/// the frame by a factor of `exp(-dt/tau)`, and scaling every leg rate by a
/// positive constant scales `p` by it too.  So a command that puts the
/// frame-open `p` at zero puts the frame-end `p` at zero, and one that makes it
/// negative keeps it negative.  What IS one frame stale is the geometry — `R`,
/// `J_v` and the ball's place on the plate are this frame's, and the rates they
/// weigh are next frame's.
///
/// **It fails safe and it can fail.**  The worst case if the constraint binds
/// every frame is a plate that holds still, which is the incumbent best.  But
/// the correction is clamped to the servo travel and retreated into the
/// holdable set like any other command, and `omega x R s` is not heave's to
/// cancel once the legs run out — a plate stranded low and tilted can still
/// carry its contact point upward through a tilt it needs for pose recovery.
/// `SimReport::contact_normal_rate` is what says whether that happened, and
/// `test_attract_mode` asserts on it rather than trusting this header.
///
/// Returns `cmd_rad` unchanged when the contact point is already level or
/// falling, when there is no lag to convert a rate into a command, and when
/// `plate.rates_trustworthy` is false — the last for the reason the contact
/// model refuses a separation on the same flag: near a singularity `J_v` is
/// arithmetic rather than physics, and inverting it would be commanding a leg
/// rate to cancel a plate motion nobody can vouch for.
///
/// `ball_s` is the ball's centre in the plate frame, the same point
/// `contactNormalRate` is defined at.
///
/// **It lives with the controller rather than with the contact model**, which
/// is why this header reaches for `PlateMotion`.  It reads the plate all the
/// way through and could sit beside `contactNormalRate` on that argument alone
/// — but what it produces is a leg COMMAND, it has to retreat that command into
/// the holdable set exactly as `legCommand` does, and it is a decision about
/// what the loop should do rather than a fact about the plant.  Putting it in
/// `ball_contact` would point the plant at the controller to borrow
/// `retreatToHoldable`, which is the dependency the wrong way round.
std::array<double, 3> holdContactDown(const TableKinematics& tk,
                                      const AutoBalanceDesign& d,
                                      const PlateMotion& plate,
                                      const std::array<double, 3>& alpha_rad,
                                      const std::array<double, 3>& cmd_rad,
                                      const Eigen::Vector3d& ball_s);

/// Advance the three rate-limited first-order servos one step, in closed form.
///
/// The servo is a first-order lag that cannot be driven faster than
/// `rate_max`, so the step is **piecewise** — ramp at the limit until the error
/// has come down to where the lag is the slower of the two, then decay exactly:
///
///     e = cmd - alpha,  and the lag takes over at |e| = rate_max * tau
///
///     |e| <= rate_max * tau :  alpha <- cmd + (alpha - cmd) exp(-dt / tau)
///     otherwise, with t_r = (|e| - rate_max * tau) / rate_max :
///         dt <= t_r :  alpha <- alpha + sign(e) rate_max dt
///         dt >  t_r :  alpha <- cmd - sign(e) rate_max tau exp(-(dt - t_r) / tau)
///
/// This header used to claim "integrated exactly" and the claim was doing real
/// work, so it is rewritten rather than quietly weakened.  **Both branches are
/// still closed form and neither can overshoot at any dt**: the ramp is capped
/// at `t_r` so it stops at the handover rather than past it, and the decay
/// approaches `cmd` without reaching it.  That matters because the plate runs
/// at a fixed 60 Hz against tau = 0.05 s — three steps per time constant, where
/// forward Euler already overshoots, and a user dragging tau below 1/30 s would
/// make Euler oscillate and then diverge.
///
/// `rate_max` comes from `TableParams::alpha_rate_max`; a non-positive or
/// non-finite one means no limit and recovers the pure exponential exactly, the
/// same way `tau <= 0` means no lag.  Why a servo has a rate at all, and where
/// 10.47 rad/s comes from, is at `kServoRateMax`.
std::array<double, 3> stepServos(const std::array<double, 3>& alpha_rad,
                                 const std::array<double, 3>& cmd_rad,
                                 double tau,
                                 double dt,
                                 double rate_max);

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
///
/// The rate limit is not a parameter here: it is `tk.params().alpha_rate_max`,
/// because it is a property of the mechanism being stepped and not a choice the
/// caller gets to make.  This is the entry point the application and every
/// harness step the plate through, so every one of them inherits it.
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
/// **It was 0.5, then 0.3, and it is 0.2 because the ball learned to bounce.**
/// Each move was a measurement rather than a taste:
///
///   - **0.5 was measured from rest at the centre of a level plate**, and the
///     demo does not open like that.  It opens tracking a circle, so the legs
///     are already displaced and the ball is already running at 75 mm/s when
///     the visitor reaches for Nudge.  A shove there lands on a plate that is
///     already working.
///   - **0.3 was that envelope measured against a ball that could not bounce.**
///     A shove hard enough to separate the ball used to end with it arriving
///     and sticking; with `kRestitution` it arrives and leaves again, and the
///     train that follows runs about `e/(1-e) = 16` times the first flight —
///     a second of a ball the plate can only reach through the horizontal
///     component of a normal impulse at each contact.
///
/// Re-measured on the same 36 x 24 grid — 36 points of the lap, 24 directions,
/// shoving a ball that is already tracking the opening circle — against the
/// contact model that ships, `holdContactDown` included, balls lost out of 864:
///
/// | speed | Nominal | Aggressive | Detuned |
/// |---|---|---|---|
/// | 0.18 | 0 | 0 | 0 |
/// | **0.20** | **0** | **0** | **0** |
/// | 0.22 | 0 | 1 | 0 |
/// | 0.24 | 0 | 0 | 0 |
/// | 0.25 | 4 | 2 | 0 |
/// | 0.30 | 10 | 35 | 0 |
///
/// The slivers thin as the shove shrinks rather than stopping at a threshold —
/// 0.22 loses one direction and 0.24 loses none, which is what a sliver looks
/// like when a 24-direction grid walks past the side of one — so the bound is
/// the last speed BELOW the first of them rather than the last clean row.
/// That is the same reading that put it at 0.3 against the old table, applied
/// to the new one.
///
/// **D11 alone would allow 0.24, and the reason it is not is #19.**  The
/// contract is that Nominal must hold every setting the sliders offer while
/// Aggressive is allowed to lose the ball — "a trap is a loss the visitor did
/// not ask for; a lesson is one they did" — and Nominal is clean through 0.24,
/// with only Aggressive's single direction at 0.22 in the way.  Read that way
/// the bound is set by the wrong row.
///
/// It stays at 0.20 because #19 made a second promise at this very constant:
/// `test_aggressive_is_no_more_fragile_than_the_shipped_tuning` asserts that
/// Aggressive survives `kMaxNudgeSpeed` from every direction, and that is what
/// makes offering the tuning behind a button honest.  Setting the bound where
/// Aggressive is known to lose a sliver would leave that assertion passing only
/// because its grid is coarser than the one the sliver was found on.  Raising
/// it to 0.24 is therefore a product decision — drop #19's comparison, or
/// exempt Aggressive from it — rather than a consequence of the table.
///
/// **What the constraint is worth, at the bound.**  Without `holdContactDown`
/// the same grid at 0.20 loses 6, 7 and 3 of 144, throws the ball up to 1.29 m
/// off the plate, and lets successive arrivals grow twelvefold; with it,
/// nothing is lost, the hop is under 110 mm and no arrival more than doubles.
/// The bound is what is left over after the constraint, not instead of it.
constexpr double kMaxNudgeSpeed = 0.20;

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
