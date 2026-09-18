# Grill brief — what the loop should do with a ball that bounces

Seed for a fresh session. Start in `/home/nds/Desktop/Merlin/04-projects/Caliburn`,
run `/grilling`, and paste everything below the line as the subject.

Use `/grill-with-docs` instead if you want the decisions to land as ADRs in
`docs/adr/` as they are made — same subject text.

---

Grill me on what the ball-balancer's control loop should do with a ball that
bounces. Read `gh issue view 23 --comments` first — the last comment has the
measurements. Branches `wip/23-round-two` and `wip/23-normal-force` in
projects/ball-balancer have the code. Don't write code this session.

## Where this came from

#23 round two asked for an elastic landing. It's built and it works: e = 0.94
(Chai et al. 2023, POM on 304 stainless), restitution on the relative normal
velocity at the contact point, sub-frame impact resolution, a termination floor
of u < g·dt/2 derived from the frame rate. Apex heights decay at e² = 0.884 on a
still plate. The physics is not what I want grilled.

What it exposed is: at e = 0.94 the bounce train runs e/(1−e) ≈ 16× the first
flight — about a second — and for that second the plate can only reach the ball
through the horizontal component of a normal impulse at each contact. There is
no rolling friction in flight and no tangential impulse at impact (a ball that
leaves while rolling without slipping keeps ω = v/r, so the contact point is at
zero slip when it lands). So control authority is duty-cycled AND small.

Consequence: 3 of the 24 path settings the sliders offer lose the ball — the
cornered shapes at their fastest laps. Circles are untouched at every size and
speed. No interface bound fixes it: losses persist at every speed cap from
0.25 m/s down to 0.10, and a 180 mm square would need a 25 s lap to be safe.

## Question 1 — what does the loop regulate on while the ball is bouncing?

Four candidates measured against a 72-direction kick sweep and a 24-setting path
sweep. All are worse than or equal to doing nothing:

1. Today: `predictedLanding` with the ball's velocity zeroed. The prediction
   collapses onto the ball at every impact and springs out on every rebound, so
   the loop chases a target oscillating at the bounce frequency. It phase-locks:
   the plate rises ~0.02 m/s into the ball at each arrival, feeding it ~5% per
   impact against the 12% e² removes. Square 180 mm at a 30 s lap: 1830 airborne
   frames of a 3640-frame run instead of settling.
2. Landing point + the ball's real velocity: 12 of 24 lost, condition 16350.
   The landing point already carries the velocity's effect — double-counted lead.
3. `predictedRest` (predict the end of the whole train; it's continuous across
   an impact, which does kill the oscillation — bounces fall 1130 → 60): 25 of
   72 lost, condition 23039. The loop chases a target a second ahead and slams.
4. Keeping the path feedforward while airborne: 5 path settings lost, not 3.

## Question 2 — should a corner put a step into the leg command at all?

Grill this one first; it may dissolve question 1. #24 made the reference
velocity deliberately discontinuous at corners ("a corner IS a step in the
reference velocity, and it is the reason the ball rounds one"). That step
reaches the actuator through K's velocity columns, the servo lag turns it into
~50 rad/s² of plate angular acceleration, and that is the single mechanism
behind every corner separation. It was harmless while the ball was glued down.

If corners stop stepping the command, the bouncing-ball control problem may only
need to hold for the nudge case — which already works: 0 lost in 72 directions,
all home, at the shipped 0.26 m/s.

## Off the table, so don't spend the session there

- Retuning e. #23 says explicitly: if the sourced value disagrees with what the
  demo needs, say so, don't quietly retune the ball. Already said, in the issue.
- The normal-force second derivative. Separately found and fixed on
  `wip/23-normal-force` — it was differencing a stepped command, overstating by
  exactly tau/dt (3× at 60 Hz, 12× at 240). That drops separations from 30-in-72
  to 2-in-720. It's validated and independently sound; it's blocked only by ~5
  near-singularity losses at the demo's edge, which is #22 territory and
  probably its own ticket.
- The hopping controller (#23's other half). It needs question 1 answered first.

## What I want out of it

A decision I can implement, or a clear statement that the bounce shouldn't land
at all and why.
