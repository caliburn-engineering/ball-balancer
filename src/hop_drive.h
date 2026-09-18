// src/hop_drive.h
#pragma once

#include "ball_contact.h"   // kRestitution, which sizes the train a throw starts

namespace caliburn {

/// The hopping controller's cycle: where the plate is in a deliberate throw.
///
/// **The hop is #23's constraint with the sign released, and nothing else.**
/// `holdContactDown` drives the contact point's normal rate `p` to zero while
/// the ball is airborne, because `u_rebound = e u + p (1 + e)` makes `p <= 0` a
/// proof that the loop cannot pump a bouncing ball.  A throw is the same
/// `heaveToContactRate` solve asked for `p > 0` instead, at a contact this
/// cycle chooses.  There is one heave actuator in the application and both
/// halves of the argument use it; see
/// [#34](https://github.com/caliburn-engineering/caliburn/issues/34).
///
/// So what is here is *when*, not *how*.  The cycle owns no kinematics, no
/// gain and no plate — it answers "what should the contact point be doing this
/// frame", and `stepSim` turns that into a leg triple through the same function
/// the constraint goes through.
enum class HopPhase {
    /// Not hopping.  The loop tracks, and an airborne ball is #23's.
    Off,

    /// Hopping is on, the ball is on the plate, and a throw is refused anyway
    /// — see `stepHop` for the two gates.  A distinct state rather than a
    /// silent `Off` because the panel has to be able to say WHY nothing is
    /// happening: a checkbox that does nothing at nine of the eighteen slider
    /// settings is a trap unless it says which nine and why.
    Waiting,

    /// All three legs going down together, buying the travel the throw spends
    /// and putting the plate below the heave it will end the hop at.
    Charge,

    /// All three legs coming up together, fast enough that the servo's own
    /// deceleration beats gravity and the plate lets the ball go.
    Throw,

    /// The ball is off the plate.  **The hop does not drive a flight** — the
    /// command goes back to `holdContactDown`, so the bounce train that follows
    /// a throw decays at `e^2` exactly as an accidental one does, and the hop
    /// height is bounded by the throw that started it rather than by a sweep.
    Flight,
};

/// How big a hop is, and how long each stroke lasts.  Both derived.
struct HopSizing {
    /// The contact point's normal rate the throw commands at the end of its
    /// stroke, in m/s.
    ///
    /// **Bounded below by the servo, not chosen.**  A ball is let go when
    /// `N/m <= 0`, and under a held command the leg's own lag decelerates it at
    /// `alpha_dot / tau`, so a plate rising at `z_dot` is accelerating downward
    /// at `z_dot / tau`.  Separation therefore needs `z_dot > g tau` — 0.49 m/s
    /// against the shipped 0.05 s lag — and a throw slower than that is not a
    /// small hop, it is a plate lifting the ball and setting it down again.
    /// Measured: at `1.1 g tau` the ball never leaves the plate at any point of
    /// a lap.
    ///
    /// The margin above that floor covers what the bare threshold leaves out —
    /// a tilted plate holds the ball by `g n_z` rather than by `g`, the
    /// rotational term under an off-centre ball carries the contact point
    /// either way, and the retreat into the holdable set gives some of the
    /// commanded rate back.  It is bounded ABOVE as well, and tightly: see
    /// `kHopRiseMargin`.
    double rise_rate = 0.0;

    /// How long each of the two strokes may last, in seconds, and so how long
    /// the commanded rate takes to ramp from nothing to `rise_rate`.
    ///
    /// **The ramp is what makes the plate carry the ball rather than hit it**,
    /// and the stroke is how long it is.  `heaveToContactRate` commands a rate
    /// by setting a position error of `tau * rate`, so a rate asked for flat
    /// arrives inside one frame: the plate steps straight past the `g tau`
    /// release threshold, the ball is let go at the velocity the plate had
    /// BEFORE the step, and the plate then runs into it.  That is an impact,
    /// `u_rebound = e u + p (1 + e)` applies, and the hop comes out at roughly
    /// `(1 + 2e)` times the rate commanded.
    ///
    /// Ramped instead, the plate accelerates under the ball and the ball rides
    /// it up until the plate can outrun gravity, at which point it simply stops
    /// being pressed.  Measured on the opening circle, against the same
    /// `rise_rate`: **one lag gives a 53 mm hop and gives the no-pumping
    /// constraint up by 0.44 m/s; two lags give 24 mm and give it up by 0.0002**
    /// — which is the constraint holding, not yielding.  Three lags measure the
    /// same as two and spend half again as much travel, so two is where it
    /// stops buying anything.
    double stroke_s = 0.0;

    /// How far the charge stroke takes the plate down, in metres — the area
    /// under the ramp, `rise_rate * stroke_s / 2`.
    ///
    /// **The workspace cost of a hop, and the whole of it.**  The throw is cut
    /// short the frame the ball leaves and the flight is `holdContactDown`'s, so
    /// this is the largest heave excursion a hop asks for.  Reported rather than
    /// left to be multiplied out, because it is the number the tilt authority
    /// has to be quoted against — see `test_hop_drive`, which quotes it.
    double drop_m = 0.0;

    /// How long a throw keeps the ball off the plate, in seconds: the first
    /// flight and the whole bounce train after it, `2 u / (g (1 - e))`.
    ///
    /// **`e = 0.94` is what makes this the governing number.**  A throw at `u`
    /// buys a first flight of `2 u / g`, and the train that follows is
    /// `1 / (1 - e)` times it — about seventeen.  So a 25 mm hop is not a
    /// quarter-second event, it is two and a half seconds of a ball the plate
    /// can only touch at seventeen instants.  That is also the answer to which
    /// contact to hop on: there is no choice to make, because the next throw
    /// cannot be armed until the train has ended.
    double flight_s = 0.0;

    /// The fastest the ball may be crossing the plate for a throw to be armed,
    /// in m/s.
    ///
    /// **A throw hands the ball a ballistic coast, and the coast has to fit on
    /// the plate.**  `room / flight_s`, where `room` is `R_table - r_ball` —
    /// the ball's own rim.  Counting the plate's steering during the train as
    /// zero is the conservative half of that: the plate does get seventeen
    /// normal impulses to aim, and each one turns the ball.  Taking the ball as
    /// starting from the CENTRE is the optimistic half, and is why this is a
    /// bound checked by a sweep rather than a proof — the position-inclusive
    /// version refuses to hop at all on the 180 mm paths, which measure
    /// perfectly safe.
    ///
    /// It comes out at 0.112 m/s on the shipped plate — 0.28 m of room against
    /// a 2.50 s train — and the sweep the derivation was checked against agrees
    /// to the resolution it was run at.  Over the 18 shape-size-lap settings
    /// the sliders offer: a gate that lets everything through loses the ball on
    /// **four**, all four of them a polygon at its fastest offered lap where
    /// the ball is crossing the plate at a quarter of a metre a second; a gate
    /// at 0.10 loses it on none; and the derived 0.112 that ships loses it on
    /// none either — `test_hop_drive` is what pins that, at the shipped value
    /// rather than at the exploratory one.
    ///
    /// What the refused settings get is the tracking they had before, bit for
    /// bit, which is the honest thing for a demo to do with a feature that does
    /// not fit: nine of the eighteen never hop.
    double max_ball_speed = 0.0;
};

/// The margin on `g tau` the throw's rate is ramped up to.
///
/// **Bounded on both sides by measurement, and the window is narrow.**  Swept
/// against a ball tracking the opening circle, over stroke lengths of one and a
/// half, two and three lags:
///
/// | margin | hop | tracking, mean | residual rise in flight |
/// |---|---|---|---|
/// | 1.1 | the ball never leaves | 3.1 mm | — |
/// | 1.4 | 21 mm | 3.7 mm | 0.0002 m/s |
/// | **1.5** | **24 mm** | **4.0 mm** | **0.0002 m/s** |
/// | 1.7 | 62 mm | 5.0 mm | 0.42 m/s |
/// | 2.2 | 112 mm | 14.1 mm | 0.69 m/s |
/// | 3.0 | the ball is thrown off the plate | — | — |
///
/// The step between 1.6 and 1.7 is the release turning from a carry into a
/// strike: past it the ramp crosses `g tau` by more than a frame's worth of
/// rate, the plate runs into the ball it has just released, and the `(1 + e)`
/// term triples the hop.  That is also where the loop starts running its servos
/// out of speed during the flight and giving the no-pumping constraint back.
/// 1.5 sits below the step with the whole of it in hand.
inline constexpr double kHopRiseMargin = 1.5;

/// How many servo lags a stroke lasts.  See `HopSizing::stroke_s` for the
/// measurement that puts it at two rather than at one or three.
inline constexpr double kHopStrokeLags = 2.0;

/// How big a hop is on this plate: gravity, the servo lag it is thrown by, and
/// the room the ball has to come down in (`R_table - r_ball`).
///
/// Every argument is the plant's, so a slower servo asks for a faster throw and
/// gets a higher hop — which is the correct way round, since a lag that takes
/// longer to stop is a lag that has to be outrun by more — and a smaller plate
/// refuses to hop a ball that is moving as fast.
HopSizing hopSizing(double gravity, double servo_tau, double room_m);

/// Where the cycle is, and what it has done.  One frame hands this to the next,
/// so it lives in `SimState` beside the path phase.
struct HopCycle {
    HopPhase phase = HopPhase::Off;

    /// Seconds spent in `phase`, so far.
    double phase_s = 0.0;

    /// Throws that actually let the ball go, since the cycle was switched on.
    /// A throw whose window expires with the ball still on the plate is not
    /// counted — it is the failure this number exists to make visible.
    int thrown = 0;

    /// The plate's heave when this hop was armed, in metres.
    ///
    /// **What makes `HopSizing::drop_m` a bound rather than an aspiration.**
    /// The strokes command a RATE, so a charge that ends on a timer ends
    /// wherever the plate happened to get to — and a throw that fails to let go
    /// then recharges from there, and the next from there again.  Measured
    /// before this field existed, the plate walked 101 mm below its home heave
    /// on a path whose stated cost was 37, and lost the ball on 4 of the 18
    /// slider settings it was swept over.  A hop is bounded to
    /// `[z_arm - drop_m, z_arm]` instead, and the plate comes back to where it
    /// started whether or not the throw worked.
    double z_arm = 0.0;
};

/// What the cycle is shown of the frame it is deciding about.
///
/// Bundled rather than passed as six positional arguments, three of them
/// adjacent bools: `stepHop(size, cycle, true, false, true, 0.0, z, dt)` is a
/// call nobody can read and any pair of which can be swapped silently.
struct HopFrame {
    /// Whether the hop is switched on at all.  False RESETS the cycle rather
    /// than pausing it, so switching the hop off mid-throw does not leave a
    /// half-charged plate to resume into.
    bool enabled = false;

    /// Whether the ball is off the plate.  The flight is #23's, not the hop's.
    bool airborne = false;

    /// Whether the TRACKING command — the one the gain produced, before any
    /// heave correction — is inside the servo travel and inside the holdable
    /// set.  See `stepHop` for why a loop on its stops may not throw.
    bool loop_has_room = false;

    /// How fast the ball is crossing the plate, in m/s, against
    /// `HopSizing::max_ball_speed`.
    double ball_speed = 0.0;

    /// The plate's heave, in metres.  The only thing the cycle knows about the
    /// mechanism, and what the strokes end on — see `HopCycle::z_arm`.
    double plate_z = 0.0;
};

/// What the cycle wants of the plate's contact point this frame.
struct HopDemand {
    /// False leaves the command alone — the loop tracks, and `holdContactDown`
    /// constrains it if the ball is airborne.  The hop never drives a flight.
    bool drive = false;

    /// The `p` to command, in m/s: negative on the charge, positive on the
    /// throw.  Meaningless while `drive` is false.
    double contact_rate = 0.0;
};

/// Advance the cycle one frame and say what it asks of the contact point.
///
/// **Once per landing, never inside a bounce train**, and that is the answer to
/// "which contact to hop on".  A throw is armed by the ball being ON the plate,
/// and the train that follows one is about seventeen times the first flight at
/// `e = 0.94` — so hopping at every impact would be a different feature, and a
/// pumped one.  The flight belongs to `holdContactDown` from the frame the ball
/// leaves to the frame it lands, which is what bounds the hop height by the
/// throw rather than by a measurement.
///
/// **Two gates on ARMING are what stop the hop throwing the ball off the
/// plate**, which is #23's own open question about this feature and the failure
/// mode #23 was opened about:
///
///   - **`HopFrame::ball_speed` against `HopSizing::max_ball_speed`.**  A throw
///     hands the ball two and a half seconds of coasting, and the coast has to
///     fit on the plate.  This is the one that does the work: it is what takes
///     the sweep from four balls lost of eighteen settings to none.
///   - **`HopFrame::loop_has_room`.**  A loop already on its stops has no
///     authority to spare and heave spends authority.  Measured on the 180 mm
///     square at its fastest offered lap, hopping regardless walks the ball out
///     to 200 mm and off the rim with `saturated` and `clipped` set on every
///     frame of the way out.
///
/// Both gate arming rather than the stroke in progress: a hop that has already
/// charged is carried through to its throw, because a plate left low is worse
/// than a plate that finished what it started.
///
/// The strokes end on TRAVEL rather than on their timers — see
/// `HopCycle::z_arm` — and the timers are the backstop for a plate that cannot
/// deliver the travel, so that a hop against the servo stops stalls rather than
/// marching.
HopDemand stepHop(const HopSizing& size, HopCycle& cycle,
                  const HopFrame& frame, double dt);

}  // namespace caliburn
