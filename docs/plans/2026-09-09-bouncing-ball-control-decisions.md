# Decisions — what the loop does with a ball that bounces

Outcome of the `/grilling` session seeded by
[`2026-09-08-grill-bouncing-ball-control.md`](2026-09-08-grill-bouncing-ball-control.md).
No code was written. Eighteen decisions, six tickets, one reversal of #24 and
one correction to a claim made inside this session.

Subject: [#23](https://github.com/caliburn-engineering/caliburn/issues/23).

---

## The finding that reframes the brief

Restitution acts on the relative normal velocity at the contact point, so for a
ball arriving at speed `u` onto a plate whose contact point has normal velocity
`p`:

```
u_rebound = e·u + p·(1 + e)
```

**The sign of `p` is the whole of the pumping question, and it is a hard result
rather than a tuning.** `p <= 0` gives `u_rebound <= e·u`, so the apex ratio is
at most `e²` — always. A rising plate always adds energy, at `1.94x` its own
speed at `e = 0.94`.

This is why all four measured candidates lost to doing nothing. They were not
four control laws. K's output is a leg *triple*, which is simultaneously tilt
(differential) and heave (common-mode); those do completely different jobs while
the ball is airborne — tilt aims the normal impulse at landing, heave sets its
magnitude — and feeding a single target through K leaves their relative phase
uncontrolled. The four candidates were four different ways of **not saying
anything about `p`**.

The brief's Question 1 is therefore not "which target" but "what constrains one
scalar".

---

## Decisions

### On the corner (the brief's Question 2)

**D1 — A corner must not put a step into the leg command.** Two fixes, in this
order, because they bind in different regimes: reference shaping is an
*acceleration* limit and the servo rate limit is a *rate* limit, and separation
at a corner is an acceleration phenomenon.

**D2 — Fillet the path geometry; do not slew the velocity signal.** Slewing
`pathVelocity` without touching `pathPoint` would make the reference velocity
stop being the derivative of the reference position — a new lie in exactly the
place we are removing one, and a direct violation of `stepPath`'s stated reason
for returning both together. A circular blend of radius `v²/a_max` keeps
`v = dp/dt` true by construction.

**D3 — `pathOutline` draws the fillet.** A drawn square with a filleted
reference is a drawing that lies about the path. `pathLength`, `minPeriod` and
`phaseNearest` all move with it.

**D4 — `a_max = (5/7)·g·sin(theta_max)`, derived from feasibility.** The most
acceleration the plate can give the ball. A reference asking for more than the
ball can produce is infeasible by construction — the same principle
`kMaxSetpointSpeed` and `kMinLapSeconds` already express.

**D5 — `theta_max` is the largest tilt whose Jacobian condition stays under 20,
computed by sweeping `tk_.condition_number`, not written as a literal.** Reuses
the threshold this repo has already argued for (`kRatesUntrustworthyAbove`,
the UI's own "Poor" line) rather than inventing a second opinion about when the
mechanism is in trouble. It also makes `a_max` a derived property of the plant,
so the leg-length sliders move it — the behaviour every other bound already has.
The workspace maximum is explicitly *not* the answer: the workspace edge is
precisely where the condition number blows up.

For scale, `(5/7)g = 7.007`, and Square 180 mm at its 4.07 s floor runs a
254.6 mm side at 250 mm/s:

| theta_max | a_max [m/s²] | fillet [mm] | % of side |
|---|---|---|---|
| 10° | 1.22 | 51 | 20% |
| 15° | 1.81 | 34 | 14% |
| 20° | 2.40 | 26 | 10% |
| 30° | 3.50 | 18 | 7% |

A fillet that large is the feature, not the cost. The square at that lap time
genuinely cannot be tracked with sharp corners; today the reference lies about
that and lets the plant discover it by throwing the ball. The fillet is
`v²/a_max`, so it shrinks to under a millimetre on a 30 s lap and grows as the
lap tightens — #24's bandwidth argument made visible *in the target* instead of
inferred from the ball's overshoot.

**D6 — The servo gets a rate limit, cited, living in `TableParams` beside the
travel limits.** ~10.5 rad/s, the class a high-torque digital servo at
0.1 s/60° actually delivers. A first-order lag has unbounded initial rate; a
real servo does not, and the plant is missing a physical property rather than
needing a controller special case.

It engages at `|cmd − alpha| > rate·tau = 0.525 rad = 30°`, so it will rarely
bind at a corner — it is insurance for saturated kick recoveries and the
Aggressive preset. The property that makes it worth having: **while
rate-saturated, `alpha_dot` is constant, so `alpha_ddot = 0` exactly** — the
`c_ddot` and `omega_dot` terms that decide separation vanish precisely when the
loop is slamming hardest. At disengage `alpha_ddot` jumps to `−alpha_dot/tau`,
the same magnitude it would have had unlimited, but later and from a smaller
error. Strictly better, never worse.

**D7 — It goes in the plant (`stepServos` / `stepServosOnPlate`), so every
harness inherits it, and K is unchanged.** Same reasoning as the travel clamp: a
saturation outside the design model, accepted and documented, not a redesign
trigger. `stepServos`'s "integrated exactly" claim becomes piecewise — ramp at
the limit until `|cmd − alpha| <= rate·tau`, then decay exactly — and the header
gets rewritten rather than quietly weakened.

### On the airborne law (the brief's Question 1)

**D8 — Constrain `p`, do not hunt for a fifth target.** While airborne, saturate
the common-mode leg motion so the contact point never rises: `p <= 0`.
`predictedLanding` stays as the horizontal target and tilt authority is
untouched. Non-pumping becomes a proof rather than a measurement.

This is the same decomposition #23 already identified for the hopping
controller — "tracking is the *differential* part of the leg triple and hopping
is the *common-mode* part" — arriving a ticket early. It fails safe: if the
constraint binds constantly, the worst case is the still plate, which the brief
establishes is the incumbent best.

**D9 — A clean path sweep is necessary and not sufficient; the phase-lock is a
defect in its own right.** The fillet stops corners *causing* separation. It
does nothing to the law that runs once a separation happens, and kicks and the
Aggressive preset still separate the ball by design. A loop that injects energy
into a bouncing ball at the bounce frequency is wrong in the same way the
differenced `omega_dot` was wrong: survivable at the shipped tuning, and still
an artefact rather than physics. #23 closes on the sweep **and** the no-pumping
test.

**D10 — The no-pumping test is the apex ratio, asserted against `e²` over the
72-direction sweep with the loop closed.** The bounce was validated exactly this
way on a still plate ("apex heights decay at `e² = 0.884 ± 2%`"); the closed-loop
test is the same measurement with the loop running, and the passing condition is
that closing the loop does not make it worse. Bounce count and airborne-frame
count are secondary assertions — they are symptoms, the apex ratio is the
mechanism, and it is the one that fails loudly if a rising plate is ever
reintroduced.

**Known risk, recorded rather than resolved.** The plate descends at separation
*because* it accelerated away from the ball, so forbidding it to rise for a ~1 s
train may strand it low and tilted, and pose recovery through tilt still moves
the contact point vertically via the `omega x Rs` term. If that conflict binds
hard, the fallback is the explicit split — decompose `u` into differential and
common-mode, regulate the differential on `predictedLanding`, hold or constrain
the common-mode for the flight. Same idea, more machinery.

### On what the demo owes a visitor

**D11 — Nominal must hold every setting the sliders offer. Aggressive is allowed
to lose the ball, because that is what the visitor asked for.** This makes the
two positions already in the tree consistent instead of contradictory: a *trap*
is a loss the visitor did not ask for (`setpoint_path.h`: "a demo whose controls
include a setting that breaks it is not offering a choice, it is offering a
trap"), a *lesson* is one they did (#19's Aggressive preset). It also gives the
corner work a falsifiable pass mark — 24-of-24 under Nominal — rather than
"fewer losses".

**D12 — The contract is dated, not retroactive.** D11 binds the ticket that
closes #23, not every commit on the way there.

**D13 — A lost ball is shown and it is recoverable.** Banner plus a
`plot_state_.markers` entry at the instant, plus the sim pauses; the existing
Auto-reset checkbox stays as the escape hatch and defaults **off**. Today
`ball_auto_reset_` defaults true and `resetBall()` silently teleports the ball to
centre with no cause given, which is #22's failure mode exactly — the demo looked
fine and was not. The plot marker is the part that makes it legible rather than
merely announced.

**D14 — One cause is named, by a precedence derived from causality and pinned by
a test.** In order: rates untrustworthy (the arithmetic cannot be believed at
all) → workspace-clipped (the mechanism refused the command) → separated (the
plate left the ball) → saturated (the servo ran out of travel) → else rolled off
(nothing failed; the ball went too far). Each is upstream of the next, so the
precedence is a claim about mechanism rather than a ranking of severity. The raw
flags go in the marker's tooltip.

### On evidence

**D15 — Re-measure the 24-setting path sweep on the analytic baseline, and commit
it as a test.** The 3-of-24 number has gated two sessions and exists only in
prose. `dccebc5` is a single commit carrying both the bounce and the analytic
normal force, so which estimator the table stands on is not recoverable from the
branch — and the demo losses and the estimator bug had the *same* trigger, every
corner of every lap. If it returns 0-of-24, the brief's Question 2 is answered by
deletion.

**D16 — Re-measure `kMaxSetpointSpeed` in the same sweep, and let it move.**
0.25 m/s was measured against a plant that slammed the legs at every corner, and
its own header's "loses it entirely beyond about 400 mm/s" is contaminated by the
delta function the fillet removes. If it comes back higher, that is the feature
restoration the 25 s lap problem was asking for. If it does not move, the header
should say what the cap is actually about.

**D17 — Collapse the five-way duplicated sim loop before writing any sweep.**
Not on tidiness grounds. The code review already left it deliberately — "fair,
it is what let finding 2 hide" — and this session then spent most of its length
unable to answer "which estimator was that measured against". We are about to add
four more harnesses that each drive the same loop the app drives. Doing it
afterwards means writing five more copies and deleting them.

**Correction to a claim made in this session.** `a_max` does *not* subsume
`kMaxSetpointSpeed` and `kMinLapSeconds`. Centripetal feasibility on the largest
circle permits `sqrt(a_max·R)` = 657 mm/s at `a_max = 2.4`, far looser than the
250 mm/s that was measured to be necessary; on a 20 mm circle `kMinLapSeconds`
binds at 63 mm/s against feasibility's 219. The principle is shared, the binding
constraint is not. All three bounds stay, and `a_max` is in practice a **corner**
bound — a corner being the only place on these paths with infinite curvature.

### On landing it

**D18 — Split the analytic normal force out and land it now.** It is
independently sound, validated against a central difference on a smooth drive,
and blocked only by ~5 near-singularity losses that are already called #22
territory. It is also the baseline every number above needs. Holding a correct
fix hostage to an unrelated failure is what entangled the two findings in one
commit and destroyed the measurement baseline in the first place. The
near-singularity losses become a **release-blocking** ticket, not a backlog item,
and the committed sweep lands *reporting* its number rather than asserting
24-of-24, tightening when that ticket closes.

---

## Work order

1. **Land the analytic normal force.** Rebase `wip/23-normal-force` clean, merge.
   File the near-singularity losses (3 path settings, 1-in-720 kicks, condition
   2 300–37 000) as a release-blocking ticket against #22.
2. **Collapse the duplicated sim loop.** One shared step driven by the app and
   every harness. No user-visible change.
3. **The reference must be feasible** — fillet, `a_max`, `theta_max` sweep, drawn
   outline, `setpoint_path.h` rewritten. New ticket; quotes #24 and states the
   reversal.
4. **Servo rate limit.** Own ticket, plant fidelity, independent criterion.
5. **The loss is shown.** Own ticket; banner, marker, pause, precedence test,
   Auto-reset defaulted off. The only one a visitor meets directly.
6. **#23 closes** on: the 24-setting sweep re-measured and committed, 24-of-24
   under Nominal, and the apex-ratio no-pumping assertion with `p <= 0` enforced.
7. **Hopping controller** — own ticket, blocked on #23.

`#24` is **not** reopened. Its deliverable shipped and works; a later finding
overturned one of its premises, which is what the new ticket records.

---

## What is reversed, and the words being reversed

Round two of #23 reversed the inelastic-landing call by quoting the reasoning it
overturned rather than quietly dropping it. Same discipline here.

**`setpoint_path.h`, `pathVelocity`, from #24:**

> Undefined for an instant at each corner, where the path's velocity is genuinely
> discontinuous; the value returned there is the edge being left. That is honest
> — a corner IS a step in the reference velocity, and it is the reason the ball
> rounds one.

The first half stands: the *path's* velocity is genuinely discontinuous. The
second half does not. The corner-rounding comes from the position error against
closed-loop bandwidth, not from the velocity step; the step's only other effect
is to hand the actuator an impulse, and a reference that demands infinite
acceleration is a modelling error rather than a demonstration. What made this
visible is that it was harmless while the ball was glued down.

---

## Deliberately not decided

- **`e` stays 0.94.** Off the table by the ticket's own rule: report the
  disagreement, do not retune the ball.
- **The hopping controller's design.** Blocked on #23. Worth recording now: if
  D8 holds, the hopping controller is *the same actuator with the sign
  released* — deliberately drive `p > 0` at a chosen contact instead of clamping
  it to zero. The constraint work is not a detour before the feature, it is the
  feature's first half.
- **`bounceFloorSpeed` against a deliberate hop.** Waits on the hopping
  controller. Met by construction today (0.34 mm at 60 Hz against a 6.6 mm
  passive hop), but the interaction is untested because nothing hops on purpose
  yet.
- **Whether `predictedLanding` moves into `ball_contact`.** Noted by the review
  as real and not a bug; unchanged by any of the above.
