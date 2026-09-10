// src/setpoint_path.h
#pragma once

#include <Eigen/Core>

namespace caliburn {

/// The shapes the ball can be asked to trace.
///
/// `Fixed` is the setpoint the balance loop has always had — a point the
/// visitor drags.  The rest move it, which asks a different question of the
/// controller: not "can you hold a ball still" but "can you make it go where
/// I say".  A circle the ball actually follows is self-evidently a working
/// controller to someone who does not read a pole-zero map, and the corner it
/// rounds off on the square is self-evidently a bandwidth limit.
///
/// See [#24](https://github.com/caliburn-engineering/caliburn/issues/24).
enum class PathShape { Fixed, Circle, Square, Triangle };

/// A closed path for the ball setpoint to run round, in the plate's own frame.
struct SetpointPath {
    PathShape shape = PathShape::Fixed;

    /// Circumradius: the distance from the centre to the corner a polygon's
    /// fillet is cut from, and for the circle simply its radius.  Sized from
    /// the corner rather than the edge so that every shape at the same setting
    /// reaches equally far out, which is what a visitor comparing them expects.
    ///
    /// A filleted polygon does not quite reach it — the blend cuts the corner
    /// off, and what is left reaches `radius_m - rho*(sec(pi/n) - 1)`.  On the
    /// 180 mm shapes that is 14 mm for the square at its fastest offered lap
    /// and a quarter of a millimetre at its slowest; for the triangle 33 mm
    /// and half a millimetre, since `sec(pi/3) - 1` is exactly one and a
    /// triangle therefore gives up its whole blend radius.  A sharper corner
    /// loses more of itself to the same blend, which is the correct way round.
    ///
    /// The corner is still where the size slider is measured from, because the
    /// corner is what the shape IS — and because sizing from anything else
    /// would stop the shapes reaching equally far at the same setting, which
    /// is the promise this field was written to keep.
    double radius_m = 0.12;

    /// Seconds for one lap.  Time rather than speed because it is what makes
    /// the shapes comparable: at the same period a triangle and a circle come
    /// round together, and the triangle is simply travelling faster to do it.
    double period_s = 10.0;

    /// The most acceleration the plate can give the ball, in m/s^2, which is
    /// what the polygons' corners are blended to be feasible against.
    ///
    /// **A property of the plant, not of what the visitor asked for.**  It is
    /// `maxBallAccel`'s answer, and `stepSim` stamps it over whatever a caller
    /// left here — so no harness can measure a reference the plate could not
    /// have followed, and none has to remember to ask.  Set it directly only
    /// when reasoning about the geometry on its own, which is what
    /// `test_setpoint_path` does.
    ///
    /// Zero is the unfilleted reference this code had before
    /// [#31](https://github.com/caliburn-engineering/caliburn/issues/31):
    /// sharp corners, and a step in the reference velocity at each one.  It
    /// is kept reachable because it is the thing the fillet is measured
    /// against, not because it is a setting anyone should ship.
    double accel_max = 0.0;
};

/// The fastest the setpoint may be driven, and the furthest out it may go.
///
/// Both measured, both against the shipped tuning with velocity feedforward.
/// The ball trails the setpoint and overshoots the corners, so the path's own
/// radius is not the radius the BALL reaches: at 180 mm and 283 mm/s the ball
/// swings to 195 mm, and the plate loses it entirely beyond about 400 mm/s —
/// a 250 mm circle at a two-second lap throws it clean off.
///
/// So the sliders are bounded rather than left to find that out.  A demo whose
/// controls include a setting that breaks it is not offering a choice, it is
/// offering a trap.
inline constexpr double kMaxSetpointSpeed = 0.25;   ///< [m/s]
inline constexpr double kMaxPathRadius = 0.18;      ///< [m]

/// And a floor on the lap time whatever the size.
///
/// A speed cap alone is the wrong bound at small radii: 250 mm/s round a 20 mm
/// circle is a half-second lap, which is an angular rate the plate cannot
/// follow however short the distance.  Two constraints because there are two
/// ways to ask for something impossible — go too far too fast, or go round too
/// often — and neither implies the other.
inline constexpr double kMinLapSeconds = 2.0;

/// The shortest lap that keeps a path under `kMaxSetpointSpeed`.
///
/// The larger of the two bounds: `length / kMaxSetpointSpeed`, which binds on
/// the big paths, and `kMinLapSeconds`, which binds on the small ones.
///
/// **Not `pathLength(p) / kMaxSetpointSpeed`, and it cannot be.**  A polygon's
/// length now depends on its lap time, because the fillet does — so asking
/// `pathLength` here would be asking about the lap the caller currently has
/// rather than the one being solved for, and the answer would move every time
/// it was applied.  The way out is that at the cap the speed is KNOWN: the
/// blend radius is `kMaxSetpointSpeed^2 / accel_max` whatever the lap turns
/// out to be, so the perimeter is, so the lap that runs that perimeter at the
/// cap is.  It is a fixed point of the circular definition, reached in closed
/// form rather than iterated to.
///
/// The fillet shortens the path, so this floor comes DOWN as `accel_max`
/// falls: the 180 mm square's fastest offered lap moves from 4.07 s to 3.85 s.
/// A path that rounds its corners has less ground to cover.
double minPeriod(const SetpointPath& p);

/// The radius the polygon's corners are blended with, in metres — `v^2/a_max`
/// at the speed the path is actually being walked at, and zero for the shapes
/// that have no corners.
///
/// **The whole of #31 is this number.**  It is the radius of the tightest
/// circle the ball can be asked to follow at this speed, so a reference built
/// out of straight lines and arcs of it is one the plate can produce: the
/// centripetal acceleration through the blend is exactly `accel_max`, and it
/// is zero everywhere else.  A sharp corner is the same expression with the
/// radius sent to zero, which is where the infinite acceleration was.
///
/// It grows as the lap tightens and shrinks as the lap opens — 33 mm on the
/// 180 mm square at its floor, 0.4 mm at a thirty-second lap — which is #24's
/// bandwidth argument made visible in the TARGET instead of inferred from the
/// ball's overshoot.
///
/// Capped at the polygon's inradius, `radius_m * cos(pi/n)`, where the blends
/// meet and the shape has become its own incircle.  Beyond that there is no
/// straight left to blend and the fillet would have to cut into the previous
/// corner.  Nothing the sliders offer comes near it — the 180 mm square's cap
/// is 127 mm against the 33 mm it asks for — so it is a definition of the
/// limit rather than a clamp anyone meets.
double filletRadius(const SetpointPath& p);

/// How far round the lap the setpoint is: 0 at the start, 1 back at the start.
///
/// **Phase, not time, and this is the whole of it.**  The setpoint used to be
/// evaluated at `t / period_s`, which reads as a pure function of the clock
/// and is not one: changing `period_s` moves that quantity by `t dT / T^2`,
/// and `t` is the entire time the simulation has been running.  Measured, at
/// 100 s, a lap change of 10.0 -> 9.5 s moved the setpoint **170 degrees** —
/// so a nudge of the slider teleported the target across the plate and the
/// loop hauled the ball after it.  The jump grew with run time and wrapped, so
/// any nudge could land anywhere.
///
/// Accumulating the phase instead makes a lap change alter only the RATE from
/// that moment on, which is what the slider says it does.  See #24.
///
/// A degenerate period does not advance rather than dividing by zero.
double advancePhase(double phase, double dt, double period_s);

/// Where the setpoint is at `phase`.
///
/// The polygons are traversed at constant SPEED, not constant angle: a corner
/// is a change of direction, and slowing into it would hide exactly the
/// behaviour the cornered shapes exist to show.  Constant through the fillets
/// too — the blend changes the direction the setpoint turns through, not the
/// rate it covers ground at.
///
/// Linear in `radius_m` for every shape at a fixed `accel_max`, which is what
/// makes the size slider safe to drag: the setpoint slides straight out along
/// the ray it was already on, same angle, bigger shape.  The fillet is the one
/// thing that is not, since it is set by speed rather than by size.
Eigen::Vector2d pathPoint(const SetpointPath& p, double phase);

/// How fast the setpoint is moving at `phase`, and where it is going.
///
/// The one place `period_s` still enters directly, because it must: the phase
/// says WHERE round the lap, and the period says how fast that is being
/// walked.  So the reference velocity DOES step when the lap slider moves,
/// while the position does not — which is correct, and is exactly what the
/// slider was asking for.
///
/// **Continuous in `phase`, which reverses what #24 wrote here.**  The words
/// being reversed:
///
/// > Undefined for an instant at each corner, where the path's velocity is
/// > genuinely discontinuous; the value returned there is the edge being
/// > left.  That is honest — a corner IS a step in the reference velocity,
/// > and it is the reason the ball rounds one.
///
/// The first sentence stands: the *polygon's* velocity is genuinely
/// discontinuous, and `accel_max = 0` still returns exactly that.  The second
/// does not.  The corner-rounding comes from the position error against the
/// closed-loop bandwidth, not from the velocity step; the step's only other
/// effect is to hand the actuator an impulse, which was harmless only while
/// the ball was glued to the plate and could not be thrown off it.  A
/// reference that demands infinite acceleration is a modelling error rather
/// than a demonstration, so the polygon the setpoint runs is filleted and this
/// is its derivative everywhere.  See #31.
Eigen::Vector2d pathVelocity(const SetpointPath& p, double phase);

/// The lap time the panel actually holds, whatever the slider asked for.
///
/// The floor moves with the size — see `minPeriod` — so this has to be applied
/// against the radius just dragged rather than the one the last frame copied.
/// The frame where those differ is exactly the frame a just-enlarged path would
/// keep the smaller path's floor and run at a speed the bound exists to forbid.
///
/// Here rather than beside the slider so that the rule has one implementation.
/// A test harness that re-derives the panel's own clamp is testing its own copy
/// of it, which is the shape of the duplication #23 was bitten by.
double clampPeriod(const SetpointPath& p, double asked_s);

/// One frame of a setpoint being driven round a path.
///
/// Zeroed by default, and that is load-bearing rather than tidiness: Eigen's
/// default constructor leaves a fixed-size vector UNINITIALISED, so a
/// `PathStep{}` standing in for "no path is driving this frame" would otherwise
/// hand the loop a reference velocity of whatever was on the stack.
struct PathStep {
    Eigen::Vector2d point{Eigen::Vector2d::Zero()};      ///< [m] this frame
    Eigen::Vector2d velocity{Eigen::Vector2d::Zero()};   ///< [m/s] this frame
    double next_phase = 0.0;                             ///< for the next one
};

/// Read the setpoint, then advance the phase — in that order, which IS the rule.
///
/// Both halves are returned together because separating them is how they come
/// to disagree: the position and the velocity handed to the loop have to be the
/// same instant, and on a polygon a phase advanced in between would straddle a
/// corner and answer with the wrong edge.  Returning `next_phase` rather than
/// mutating keeps the caller's phase the frame-open one for as long as the
/// frame needs it.
///
/// `dt` is the frame the setpoint is about to be held for, so a `Fixed` path
/// still returns a zero point and a phase that does not move.
PathStep stepPath(const SetpointPath& p, double phase, double dt);

/// The phase at which `p` passes closest to `target`.
///
/// For changing SHAPE without moving the setpoint.  Phase is not comparable
/// across shapes: the circle's phase zero is at +x and a polygon's first corner
/// is at the top, so the same phase is a quarter of a lap apart — measured,
/// **170 mm** between a 120 mm circle and the square at every phase in the lap.
/// Carrying the phase across a shape change therefore teleports the target to
/// the far side of the path and the loop hauls the ball after it, hard enough
/// to drive the legs into the workspace clip.
///
/// Re-seeding from the point the setpoint is already at moves it by at most the
/// distance between the two shapes themselves — 35 mm for that circle and
/// square, and zero wherever they touch.  It is the same promise the size
/// slider makes: the target does not jump, it slides to the new path.
///
/// Exact rather than searched.  The nearest point on a circle is the radial
/// projection; on a filleted polygon it is the nearest of its straights'
/// clamped projections and its blends' clamped angular ones, and there are at
/// most four of each.
///
/// The exact centre is equidistant from the whole path and so has no nearest
/// point: it answers phase zero, for every shape, because the alternative is
/// the polygon walk picking whichever edge rounding made shortest — arbitrary,
/// and free to differ between builds.  It is one press away ("Centre setpoint",
/// then choose a shape).  A point merely NEAR the centre is still decided by
/// the arithmetic, which is correct: it genuinely does have a nearest point.
double phaseNearest(const SetpointPath& p, const Eigen::Vector2d& target);

/// The path's total length, for drawing it and for reasoning about speed.
///
/// A function of the lap time as well as the size, on a polygon: the corners
/// are blended, and the blend cuts the corner shorter than the two tangent
/// lengths it replaces.  `n*side - rho*(2n*tan(pi/n) - 2*pi)`, which is 57 mm
/// off the 180 mm square's metre at its floor and 0.7 mm off it at a
/// thirty-second lap.
double pathLength(const SetpointPath& p);

/// Points around one lap, for drawing.  A circle gets `samples` of them; a
/// polygon gets its straights exactly — a straight drawn as a chord of samples
/// is a straight drawn wrong — and its blends sampled, since they are arcs and
/// an arc has no exact polyline.
///
/// **The drawing has to carry the fillet.**  A square outlined with sharp
/// corners over a setpoint that rounds them is a picture of a path the ball is
/// not being sent round, and the visitor would read the gap as the controller
/// failing at the corner rather than as the corner not being there.
void pathOutline(const SetpointPath& p, int samples,
                 Eigen::Matrix<double, 2, Eigen::Dynamic>* out);

}  // namespace caliburn
