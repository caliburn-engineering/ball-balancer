# Ball-Balancer — Context

One application over two halves: a 3-RRS table balancing a ball under real
rolling dynamics, and the analyzer that models a plant, closes a loop around it,
and shows the consequences across Bode, Nyquist, pole-zero and step-response
views.  The repository was named `linear-analyzer` until
[#20](https://github.com/caliburn-engineering/caliburn/issues/20): it is the
survivor of the merge and carries the whole history, so the name moved to match
the product rather than the other way round.  The pre-merge ball-balancer
repository is kept on disk, unmodified, as `ball-balancer-legacy`.

## Glossary

Terms below are the project's ubiquitous language. Where a term has a plausible
synonym that is *wrong* here, the synonym is recorded with the reason — drifting
to it silently reintroduces a resolved confusion.

### Loop

A pairing of one plant output to one plant input, carrying its own gains. The
controller is a **list** of loops (`std::vector<Loop>`), not a per-channel
diagonal block — this is what lets square, SIMO and MISO plants share one data
model. A square plant with `m` loops reproduces the diagonal case.

Fan-out is legal in both directions; many-to-one is the primary MISO structure.
An empty loop list is equivalent to `ControllerType::None`.

Decided in [#2](https://github.com/caliburn-engineering/caliburn/issues/2).

### RGA (Relative Gain Array)

`Λ = G ∘ (G⁻¹)ᵀ` on the **complex** `G(jω)`, evaluated at the intended
closed-loop bandwidth — not at DC. A *pairing* diagnostic: it grades whether the
loops are matched to the right channels.

Defined only on the **square sub-matrix** of the paired channels, and therefore
only when the loop list's distinct paired outputs and distinct paired inputs are
equal and at least 2. Outside that, see **Channel Share**.

Severity comes from the RGA number `‖Λ − Π‖ₛᵤₘ`, never from `|λᵢᵢ|`.

Decided in [#4](https://github.com/caliburn-engineering/caliburn/issues/4).

### Channel Share

The fraction of one input's squared effect landing on each output, at a given ω:

```
λᵢ = |gᵢ(jω)|² / Σₖ |gₖ(jω)|²
```

Real, non-negative, sums to exactly 1. Shown wherever the RGA is vacuous or
undefined — a single-input plant, a multi-input plant with only one loop, or a
MISO fan-out whose paired sub-matrix is non-square.

**Channel Share is not a pairing measure and must never be presented as one.**
It answers "which outputs does this input reach", not "are these the right
pairings". It is also **not scale-invariant**, unlike the square RGA: it is
reported against user-supplied per-output scales (maximum allowed deviation),
and the readout states when those scales are still at their defaults.

Because it is normalised it cannot distinguish "reaches both outputs strongly"
from "reaches neither", so it is always displayed alongside `|gᵢ(jω)|` in dB.

Decided in [#8](https://github.com/caliburn-engineering/caliburn/issues/8).

> **Not "reach", "reachability", or "control authority".**
> *Reachability* is the textbook synonym for controllability, which this app
> already computes (`checkControllability`) and displays in the System Properties
> panel — the collision would be direct. *Control authority* asserts the input is
> commanded, which is false for the Quarter-Car preset, whose single input is the
> road profile disturbance.

### Structurally dead channel

A plant channel whose `|gᵢ(jω)|` is below tolerance at **every** ω on the Bode
grid — a genuine structural disconnect, as opposed to a numerical dip at one
frequency. Positive output rescaling cannot turn a zero into a nonzero at any
frequency, so this test is unaffected by the scale-dependence of Channel Share;
a threshold on the share itself, or a ratio against the column maximum, would
not be.

A loop pairing an output to an input that cannot affect it is the Ball-Balancer
failure case: `B(2,1)` / `B(3,0)` mean input 0 drives output 1's chain, so the
seeded identity pairing is dead.

Dead channels are suppressed with a stated reason rather than drawn with
misleading fallbacks — `arg(0)` renders a flat 0° phase, and the pole-zero pencil
returns `eig(A)` as coincident zeros, which reads as full cancellation.

Decided in [#8](https://github.com/caliburn-engineering/caliburn/issues/8) and
[#9](https://github.com/caliburn-engineering/caliburn/issues/9).

### Pairing grid

The controller section of the model panel: a `p × m` grid of cells, one per
plant channel, that is **simultaneously the pairing editor and the diagnostic
readout**. Clicking a cell creates or destroys the loop on that channel; each
cell shows the RGA `λᵢⱼ` (or the Channel Share, where the RGA is undefined) over
`|gᵢⱼ|` in dB. A leading column carries the per-output scale.

It is the only editor, so the loop list has no add, remove or reorder controls,
and grid order *is* loop order. One consequence: an exact duplicate pairing is
unreachable through the UI, so **#2's duplicate warning needs no surface** —
though `buildLoopController` still sums duplicates, which is correct for a loop
list built programmatically.

**Paired-ness is the cell border; severity is the cell fill.** They are separate
visual channels because a single one cannot carry both, and severity applies to
RGA cells only — Channel Share is never coloured (see above).

Decided in [#6](https://github.com/caliburn-engineering/caliburn/issues/6).

> **Not the "coupling matrix" or the "RGA table".**
> The grid outlives the RGA: on a single-input plant the same widget shows
> Channel Share, and on a 1×1 plant it is one cell with no diagnostic at all.
> Naming it after one of its readouts hides that it is first an editor.

### Plate view

The ball-balancer half of the application: the 3-RRS table, the ball on it, the
3D scene, and the Plate Control / Plate Plots panels.  It owns no window, no
ImGui context and no main loop — the analyzer's entry point owns all three, and
reconciling that was the whole of the merge.

Drawn into the **central dock node**, which is deliberately left empty so the
passthru dockspace shows the scene through it.  A window docked there would
paint over the plate.

Both targets since [#15](https://github.com/caliburn-engineering/caliburn/issues/15).
The renderer draws through the GL subset shared by OpenGL 3.3 core and WebGL2 —
vertex arrays, buffer objects, array draws — so the port was the shaders and
nothing else: `#version 330 core` and `#version 300 es` plus a precision
qualifier, supplied as a separate source string ahead of one shared body.

`glLineWidth` above 1.0 is honoured or ignored at the driver's discretion under
WebGL2, so line weight is decorative here and carries no information: every
element of the scene is told apart by colour.

### Status badges, and why nothing gets its own line

**Text that comes and goes must never change the layout's height.**

A warning drawn on its own line pushes every control below it down while it is
showing and lets them snap back when it stops.  Several of these conditions are
recomputed every frame — `balance_clipped_` and `balance_saturated_` both are —
so the panel below them does not settle into a new position, it *vibrates*, and
a slider the visitor is reaching for moves out from under the cursor.  Two
independent warnings make it worse: the controls below have three different
resting positions depending on which are showing.

Two rules, and both are needed:

- **A conditional message rides on a line that is always drawn**, appended with
  `statusBadge` (`panels/panel_utils.h`).  It costs no vertical space, it
  cannot reflow anything, and the line it attaches to is usually the reading it
  qualifies: `error: 178.2 mm  clipped  saturated` says in one line both that
  the loop is a long way off and why it is not closing the gap.  The label is a
  word or two — a badge that wraps has reintroduced the problem — and the
  sentence explaining what the condition *means* goes in the hover tooltip,
  where there is room for it.
- **Where a message is one of several alternatives, draw every branch.**  A
  trailing `else` with a disabled "in contact" or "drag to command the legs" is
  not clutter; it is what holds the line open so the panel stays still.  This
  half cannot be enforced by a helper, only by remembering it.

The block that reports the ball follows both: three lines, always — position,
velocity, and exactly one contact line out of four possibilities.  It used to
collapse to a single line when the ball went off the edge, moving everything
below it up by two rows at the precise moment the visitor was reaching for
Reset Ball.

The plot panels get the same treatment through `drawSuppressedTraces`, and
there the stakes are higher than a shifted control.  The "why is this trace
missing" line has to be drawn *above* the plot — a plot consumes
`GetContentRegionAvail()`, so anything after it is clipped out of the panel
entirely — which means a line that comes and goes does not merely move the
widgets below it.  It resizes the plot: the axes rescale and the whole curve
moves, for a reason the reader cannot see.  So the line is always there, reading
`suppressed: none` when nothing is.  The trace's NAME stays visible in the
warning colour, because a suppressed trace still has to be distinguishable from
a broken one; only the sentence saying why moves to the hover.

Not every conditional line is a defect.  A section that appears because the
visitor changed a mode — the path sliders under a trajectory, the setpoint
controls under an engaged loop — is a layout change they asked for and expect.
The rule is about text that appears because the *plant* did something.

### Tilt convention (phi/theta vs alpha/beta)

The two halves name the plate's tilt differently, and the mapping is neither
identity nor a simple swap:

```
alpha = -phi        beta = theta
```

`TableKinematics` uses `R = Ry(theta) * Rx(phi)`, so the plate's local +Y edge
*rises* with phi while its local +X edge *falls* with theta.
`RollingBallDynamics` accelerates the ball by `+g*sin(alpha)` along y and
`+g*sin(beta)` along x.  Theta therefore carries straight through and phi has to
flip.

A wrong mapping compiles, and still moves the ball — it just rolls uphill.  That
is why the convention is pinned by `tests/test_ball_sim.cpp` rather than left to
the eye, and why it lives in one function, `ballTiltFromPose`.

### Plate frame

The ball's state `[x, y, vx, vy]` is expressed in the **plate's own frame**, not
the world.  Rendering lifts it into the world with the same rotation that places
the leg attachment points: `p = c + R * (x, y, r_ball)`.

The plate is a **disc**, so the ball leaves it at a radius of
`R_table - r_ball`.  `RollingBallDynamics::on_plate` tests the inscribed
*square* and is wrong here; `ballOnPlate` is the test that matches the geometry
being drawn.

### Balance loop

The closed loop that the LQR design surface exists to produce: `u = -K(x - x_ref)`
evaluated once per frame against the **nonlinear** plate, in `auto_balance.h`.
`x` is the cascade plant's state read straight off the simulator — three leg
deviations from the linearisation point, then the ball's `[x, y, x', y']` in the
plate frame — and `u` is a leg-angle command, so nothing adapts tilt to legs
anywhere.

Engaged from Plate Control, and only while the model panel is offering a gain
that was solved against **this** plate: the right controller type, a successful
solve, a `3 x 7` gain, a ball actually being simulated, and a plant whose
geometry *and gravity* match the simulated one.  Gravity is a cascade parameter
and the plate's is fixed, so a gain designed on the moon is exactly as wrong as
one designed for longer legs and is refused the same way.  Losing any of those
drops the loop rather than freezing the last command, because a stale gain is
not a controller — and `setDesign` is the single place that drop happens, so a
draw pass cannot disagree with the step about who is driving.

The plant a cascade parameter list describes is built by `cascadeMechanism` /
`cascadeGravity` / `cascadeServoTau` / `cascadeHomeLegAngle`, and the model
builder, the application and the tests all read those same functions.  A second
transcription is how the check comes to reject two identical plates — it did,
once, over a `float` slider widening to `0.3000000119`.

Decided in [#16](https://github.com/caliburn-engineering/caliburn/issues/16).

> **The sliders and the loop cannot both be live.**
> Plate Control's three servo sliders and `u = -Kx` write the same array, so
> while the loop is engaged the sliders go read-only and keep displaying what
> the controller is asking for — the most direct view of the loop working.
> Steering moves to a **ball-position setpoint**, which is a quantity the
> design actually regulates.  A ball at rest anywhere on a flat plate is an
> equilibrium of this plant, so that setpoint enters as a reference state and
> needs no feedforward at all: the regulator is already a tracker.

That last sentence is true of a setpoint being **held**, and #24 made it worth
saying so.  A moving setpoint is not an equilibrium, and `[.., x_sp, y_sp, 0, 0]`
then claims the ball should be at the setpoint *and stationary* — so the
reference state grew a velocity.  `legCommand` takes a `BallReference`: where
the ball should be, and how fast the setpoint is travelling.  Zero velocity is
the held case and is what every caller but the trajectory driver passes, which
is why the loop went so long without one.

### Preset tunings

Three named LQR tunings — **Detuned**, **Nominal**, **Aggressive** — in the
model panel's LQR section, differing only in the ball's own state weights.  A
visitor who flips between the outer two learns more about what controller
design does in five seconds than the weighting matrices could tell them in five
minutes, which is the whole reason they exist.

> **Not the *plant* presets.**  "Preset" already means a built-in plant model
> everywhere else here — `state.preset_index`, the `Preset` combo at the top of
> the same panel, "the Quarter-Car preset".  These are **tunings**: a choice of
> Q and R against one plant, not a choice of plant.  The surface says "Tunings"
> for that reason, and the code says `LqrPreset` rather than `Preset`.  Reading
> one for the other turns "the preset changed" into two different bug reports.

| | Q on ball position / velocity | settles from 60 / -40 mm | what it shows |
|---|---|---|---|
| Detuned | 20 / 2 | 11.4 s | sluggish; never goes near a travel limit |
| Nominal | 100 / 10 | 1.8 s | the shipped default, unchanged |
| Aggressive | 150 / 2 | 0.62 s | saturates the servos the whole way in |

`R` is 1 in all three and the leg weights stay at 1, so exactly one thing
varies.  Pricing the legs *higher* buys sluggishness the position weight
already buys more legibly; pricing them lower buys a tuning that loses the ball
— `R = 0.1` against the nominal `Q` throws it off on 8 of 72 kick directions.

**`defaultLqrStateWeights` is defined as the Nominal preset**, not the other way
round, so "Nominal is what the application opens on" is arithmetic rather than
two sets of literals that have to keep agreeing.

**They are an offer, not a mode.** The weights stay editable after a preset is
chosen, and no preset is stored — `activePreset` asks the *weights* which one
they are, every frame, and answers -1 the moment a slider moves.  A stored
selection would go stale on the first drag and keep claiming a tuning the plant
no longer has.

> **"Aggressive" saturates; it does not overshoot.**
> The ticket asked for "fast but overshooting", and that is not what this plant
> does.  The servo travel limits bite before the ball ever reaches the
> setpoint, so saturation *limits the approach* rather than producing an
> overshoot — measured position overshoot past centre is under a millimetre at
> every tuning that keeps the ball.  The surface says "saturates", and the
> `saturated` and `clipped` badges in Plate Control are the honest thing to
> point at.

On the cascade the presets replace the **Reset weights** button that #16
shipped: Nominal writes exactly what it wrote, and two controls doing one thing
only raise the question of how they differ.  Off the cascade there are no
presets to reset *to*, and the button is still there.

What bounds the aggressive end is the **ball**, not the solver.  Since
[#23](https://github.com/caliburn-engineering/caliburn/issues/23) the ball can
leave the plate vertically, and a tuning that slams the legs moves the plate
out from under it rather than tilting under it — see "Contact, and the two
phases" for what an over-aggressive gain does to a ball it can no longer hold.
`Q` on position at 300 settles faster still and then throws the ball clean off
on a 0.43 m/s nudge.

150 is the fastest tuning measured that is **no more fragile than the tuning
the demo already ships**: at `kMaxNudgeSpeed`, the top of the Nudge slider,
neither loses the ball in any of 24 directions.  Both do above it, which is a
property of the plate rather than of the preset and is why the slider stops
where it does.  Pinned by
`test_aggressive_is_no_more_fragile_than_the_shipped_tuning`, so raising that
ceiling fails a test rather than quietly invalidating this paragraph.

Detuned interacts with attract mode harmlessly: it settles in 11.4 s against a
4 s kick period, so kicks would pile up — but choosing a preset is a click, and
attract mode stands down on the first click.

Decided in [#19](https://github.com/caliburn-engineering/caliburn/issues/19).

### Leg command vs leg angle

Two arrays, since the loop was closed.  `alpha_cmd_deg_` is what the sliders,
the animation or the controller **ask for**; `alpha_deg_` is where the legs
actually **are**, one first-order lag behind at the plant model's own `tau`.

The lag applies to **every** writer, not just the loop.  It is a property of
the plate rather than of the controller, and making it conditional would mean
the manual plate and the modelled plate are two different machines.  Home/Low/
High therefore command rather than teleport, and the Animation amplitude is the
amplitude *asked for* — attenuated about 1% at its default 0.5 Hz, which is
what a real servo does.  `snapServos` is the one path that writes both arrays,
because a reset is not driving anything.

Before this the slider *was* the leg angle, and closing `u = -Kx` around that
is an algebraic loop on the leg states — the gain would have been reacting to
its own command one frame earlier.  The servo dynamics the cascade model claims
had to become real in the simulator for the design to mean anything.

The lag is integrated exactly, `alpha <- cmd + (alpha - cmd)*exp(-dt/tau)`, not
by forward Euler: the plate runs at a fixed 60 Hz and tau defaults to 0.05 s,
three steps per time constant, and the tau slider goes lower still.

### Assembly mode

Which of the two solutions of the 3-RRS constraint equations a pose is.  The
mechanism is **built upward** — table above its knees — and there is a second,
**folded** root with the table lying flat on the base.  For this plate the
folded one is not an approximation or a numerical artefact: `R_ground ==
R_table` with `L1 == L2` makes `(phi, theta, z_c) = (0, 0, 0)` satisfy every
leg's length constraint exactly, at every servo angle.  It is a single pose,
the same one whatever the legs are doing.

A residual norm cannot tell the two apart, so **a converged forward-kinematics
solve is not by itself an answer about a mechanism**.  `forward_kinematics` is
a Newton solver and will return whichever root its seed is nearest;
`solve_pose` is the one that stays on the built assembly, re-solving from the
analytic level pose when the warm start lands on the folded one.

This is why issue #22 was invisible for so long.  A hard manoeuvre carried the
pose down through the knee plane, Newton settled on the folded root, and every
frame after that was seeded from a pose that was *already* a root — so it
converged in one iteration with a residual of 1e-11 and never left.  The plate
then had no tilt authority at all and every ball put on it rolled off.  Sixty
seconds in, and permanent.

The floor separating them is a tenth of the mean knee height.  Measured over
300k samples against assemblies verified by IK round-trip, the built population
runs from 0.163 to 2.00 mean knee heights and the folded root sits at zero, so
the floor clears the lowest real pose by 1.6x and the folded one by everything.
Both ends are pinned by `tests/test_assembly_mode.cpp`.

Decided in [#22](https://github.com/caliburn-engineering/caliburn/issues/22).

> **Not "the FK failed" and not "a singularity".**
> Nothing failed: the solve converged, to a real root, of the right equations.
> And a singularity is a place where the Jacobian loses rank — this is a
> perfectly well-conditioned second solution. Calling it either invites the fix
> that was already there and did not work: a tighter tolerance, or a retry
> keyed on convergence.

### Workspace vs. servo box

The three servo travel limits form a **box**; the set of leg triples the plate
can actually be assembled at is not one.  Barely half the box has any assembly,
and a *single* leg taken to either of its shipped stops with the other two at
home already has none — the mechanism's real range about the home pose is far
narrower than its servos'.

So clamping each leg to its own travel is not enough.  A per-leg clamp can hand
back a triple the plate cannot make, and a plate that cannot make its command
has no pose: it freezes at whatever tilt it last held, and a frozen tilted plate
rolls the ball off.  That was the other half of #22 — nearly two thousand frames
of a ten-minute run.

**The scatter was the folded assembly, not the controller.**  Before the fix
the loss counts jumped about with no pattern — 0.22 m/s of kick lost the ball
where 0.25 did not, and servo travel of ±35 degrees lost 61 balls in five
minutes where ±45 lost none and ±60 lost 65.  That is not something a
stabilising loop does, and it is what said the cause was discrete rather than
dynamic: whether a particular trajectory happened to carry the pose down
through the knee plane on some frame is a yes-or-no event, and everything
downstream of a yes was already broken.  Travel limits and kick magnitude
changed *which* trajectories crossed, not how much authority the loop had.

With the assembly pinned, the response is monotonic in the disturbance the way
it should always have been: at 360 swept directions the shipped tuning loses
none at 0.35 m/s, 2 at 0.42, 14 at 0.50 and 259 at 0.60.  A curve, not a
scatter — and one that says something true about the plant.

`retreatToWorkspace` is the one answer, used at both seams: pull the target back
toward a triple known to be assemblable until it is too.  Giving up magnitude
and keeping direction is what saturation ought to do, and the returned point is
always one that was actually tested, so the workspace does not have to be convex
for it to be right.

Both seams are needed, which is the part worth remembering.  `legCommand` clips
the command against the level pose; `stepServosOnPlate` clips the *step*
against where the legs already are.  Clipping only the command leaves the
straight line the first-order lag travels along, and that line can leave the
workspace even when both of its ends are inside it.  Two kick directions in
every 360 were exactly that: one frame, mid-flight, with no assembly.

The clip is not a rare corner.  Under the shipped tuning a single attract kick
saturates the legs from every direction tested and is workspace-clipped from
177 of 180 of them; under a deliberately aggressive tuning, from all 180.  It
runs on essentially every disturbed frame.

Decided in [#22](https://github.com/caliburn-engineering/caliburn/issues/22).

> **The plate has less authority than the linear design believes.**
> The cascade model knows about servo travel and nothing about the workspace,
> so the gain will ask for tilts the mechanism cannot make, and against a
> disturbance large enough it will lose the ball rather than recover it.
>
> There is no clean ceiling: the slivers of direction where it cannot come
> back get narrower as the disturbance shrinks rather than stopping at a
> threshold.  At 5760 swept directions, 0.29 m/s of ball velocity loses six,
> 0.28 loses two, and 0.26 loses none.  Attract mode's kick sits at 0.26 for
> that reason, and nobody has proved there is not a narrower sliver still.
> Designing a gain that respects the workspace is not done, and is the thing
> that would replace this argument with a guarantee.

### Attract mode

The application driving itself for a visitor who has not arrived yet: the
balance loop engaged as soon as a gain exists, and the ball already tracing a
120 mm circle at a ten-second lap.  Named after the arcade cabinet's demo reel,
which solves exactly this problem.

A motionless canvas is indistinguishable from a broken build, and a still
balanced ball reads as a photograph of one.  A visitor gives the page about ten
seconds and will not go looking for a play button.

**This used to be a disturbance schedule, and the change is worth recording
because the original reasoning was sound.**  The demo opened with the ball 72 mm
off centre and kicked it every four seconds, on the argument that what
distinguishes a control system from an animation is *recovery*: something
disturbs the ball and the loop puts it back.  In practice it read as a fault.

Two reasons, and the first is the one that matters:

- **A kick is legible only against a still baseline.**  A ball parked at the
  centre, nudged, returning — that is a story.  But the still moment between
  kicks is also the moment the page looks like a static image, so the demo
  spent most of its time looking broken in order to make the other part legible.
  A circle needs no baseline: at every instant the ball is somewhere it was
  told to be and the plate is visibly working to keep it there.  The loop is
  not demonstrating that it *can* respond, it is demonstrating that it *is*
  responding, continuously — which is also the more honest claim, since it is
  the one running every frame.
- **The opening displacement was a step input.**  Handing a state feedback
  72 mm of error at `t = 0`, with the legs exactly at home and no servo history
  to smear it, is a step: measured, the legs swung **20.7 degrees apart within
  three frames** and the table dropped 4.5 mm before settling inside 250 ms.
  Correct, and it looked like the mechanism glitching.  Starting the ball ON
  the path at the path's own velocity makes both the position and the velocity
  error zero at `t = 0`; the same measurement then gives a **3.0 degree**
  opening swing with the table height not moving at all.  Pinned by
  `test_the_opening_does_not_slam_the_legs`.

What is given up is stated plainly, because it was a real acceptance criterion:
**the demo no longer shows disturbance rejection unprompted.**  The visitor can
still see it — the Nudge buttons are right there — and the property is still
pinned by `tests/test_attract_mode.cpp`, which sweeps a 0.26 m/s disturbance
over all 720 directions and checks the ball comes home from every one.  It is
no longer *performed*.  That sweep also moved: the kick used to be production
code the demo delivered, and is now a test fixture, because what it tests is a
property of the plant and the gain rather than of the demo.

Two decisions carry the opening:

- **The opening state is a state, not a performance.**  It is set once in the
  constructor and then the visitor owns it.  The old mode had to watch for a
  visitor arriving so it could stand its kicks down, and never resume; a circle
  has nothing to stand down, so the watching went with it.
- **The circle rather than a cornered shape.**  Its tracking error is smooth,
  so the opening reads as competence rather than as the ball stumbling at every
  corner.  The corners are the better demonstration and are one dropdown away —
  but they are an argument the visitor should choose to hear, not the first
  thing they see.  120 mm at a ten-second lap is 75 mm/s, a third of the speed
  cap, tracked to 3.7 mm.

`setDesign` is where the loop first closes, and it has to be: on the opening
frame the LQR solve has not run, so a `balance_engaged_` set true in the
constructor would be cleared by the stale-gain drop and never set again.
Engaging on the first *usable* design is what "already stabilising at load"
actually amounts to.  It fires **once** — without the latch it would re-engage
a loop the visitor had deliberately dropped, on the very next frame, which is
the demo arguing with the person using it.

The application also opens with `ControllerType::LQR` selected and the
closed-loop trace visible, since a page whose whole claim is that the loop is
closed should not open its pole-zero map on the open-loop poles.

Decided in [#17](https://github.com/caliburn-engineering/caliburn/issues/17),
revised for the opening path in
[#24](https://github.com/caliburn-engineering/caliburn/issues/24).

> **Not "demo mode", "idle mode" or "screensaver".**
> All three name a state the application is *stuck in* until something releases
> it, and two of them imply a substitute for the real thing — a canned loop
> playing where the product would be.  Nothing here is canned: the plate, the
> gain and the rolling ball are the same ones the visitor gets, and the only
> difference is who is choosing the setpoint.  "Attract" also says what the
> mode is *for*, which is the test any addition to it has to pass.
>
> Prose may still call the running page "the demo" — that is the artefact, not
> the mode.

### Contact, and the two phases

The ball is no longer glued to the plate.  `RollingBallDynamics` still models
the rolling half — it is a good model of a ball in contact — but whether the
ball IS in contact is a question it could not ask, because it had no normal
force and no vertical state.

`ball_contact.h` asks it.  `N/m = g(n·z) + n·c̈ + n·[ω̇×(Rs) + ω×(ω×(Rs))] +
2n·(ω×Rṡ)` — gravity's share along the normal, the plate being driven up or
down under the ball, the plate's rotation swinging the contact point, and
Coriolis.  **`N/m ≤ 0` is separation**: a surface can push a ball and never
pull it.

**The rates in that expression are analytic, and most of the hopping this
section used to report was the estimator rather than the plate.**  `c̈` and `ω̇`
were a finite difference of `ċ` and `ω` across consecutive frames — a defensible
estimator for a smooth signal, and these are not smooth.  The leg rate is
`(cmd − α)/τ`, and `cmd` is a zero-order hold: it *steps* every time the loop
changes its mind, and differencing a step gives `1/dt`.  The overstatement is
exactly `τ/dt` on a stepping frame — 3× at 60 Hz and 12× at 240 — so it grew as
the timestep shrank rather than converging.

The servo lag differentiates in closed form instead.  `cmd` is held across the
whole frame, so `α̈ = −α̇/τ` exactly, and the same Jacobian that carries `α̇` to
the pose rates carries it to the pose accelerations.  What is left to difference
is how `J_v` and the tilt-rate map are themselves changing, and those are
functions of the leg angles and the pose: continuous, neither of them stepping.
Both history terms are load-bearing — dropping them gets `ω̇` wrong by 13% and
`c̈` wrong in *sign*.

Measured at the corner of a square path on a thirty-second lap, a setpoint
crawling at 24 mm/s and about as gentle as this demo gets:

| | `\|ω̇\|` | `N/m` | outcome |
|---|---|---|---|
| differenced | 101 rad/s² | −7.3 | ball launched at 0.26 m/s |
| analytic | 50 rad/s² | −0.6 | contact held |

A corner is a step in the reference velocity by construction (`pathVelocity`
says so), so the artefact fired on **every corner of every lap**.

So separation is real, and it is *rare*.  Under the shipped tuning at a 0.26 m/s
disturbance the ball separates in **2 of 720 directions**, peaking at **2.5 mm**
— against the 30-in-72 and 6.6 mm this section reported while the estimator was
differenced.  The plate does heave hard when the legs do (`z_c = 0.30·sin α`, so
a 20° leg swing is 74 mm of table in a tenth of a second), and a loop rejecting a
disturbance does slam the legs.  What changed is that the plate is no longer
credited with accelerations its servos never produced.

**Distinguish the tuning from the demo here, since they parted company.**  The
shipped *tuning* can still hop, as above — but at 2 directions in 720 a visitor
pressing Nudge will almost never see it, and the claim that they would was
substantially an artefact of the differenced estimator.  The shipped *demo* no
longer does at all, because it no longer kicks the ball: tracing a gentle circle
never separates it, and `test_ten_minutes_unattended_never_loses_the_ball`
asserts exactly that — zero airborne frames in ten minutes.  Both facts are
pinned, and neither implies the other.  See *Attract mode*.

**Two phases, two frames, and the frame is not a detail.**  Rolling is natural
in the plate frame — that is where `RollingBallDynamics` integrates and where
the cascade plant's state vector is written, and neither should change because
the ball can leave.  Flight is natural in the world frame, because "only
gravity acts" is a statement about an inertial frame.  Each phase keeps its own
truth and `plateFrame` converts on demand.

A launch takes the ball's **full** world velocity, including the plate's motion
at the contact point — which is most of it when the legs are slamming.  Landing
is inelastic: a real ball bounces, but a restitution coefficient is a number
nobody here has measured, and the behaviour this exists to show does not depend
on it.

**Rolling resistance is scaled by that same normal force, not by `g`.**  It is a
normal-force effect — the contact patch deforms in proportion to how hard the
surface is pressed — so the retarding acceleration is `c_rr·(N/m)`, and
`c_rr·g` was only ever the level-plate-at-rest special case.  `derivatives` and
`stepBall` take it as an argument; `quasiStaticNormalAccel` supplies `g(n·z)`
for the one caller that cannot have the real thing, the frame where the rates
are not trustworthy, and for linearising about the flat equilibrium.

The correction is one of principle rather than of magnitude, and it is worth
being plain about that.  It moves the separation count not at all — the same
directions either way, and the peak hop by 0.0001 mm — because separation is
decided by `normalAccel`, which friction does not enter, and because near
separation `N` is small and so is anything scaled by it.  What it does change
is the dead band's derivation, from `asin(c_rr)` to `atan(c_rr)`; see below.

Decided in [#23](https://github.com/caliburn-engineering/caliburn/issues/23).

> **Three suites are red on this tree, deliberately, and all three fail under
> the *Nominal* tuning.**
> `test_auto_balance`, `test_trajectory` and `test_attract_mode`.  The analytic
> normal force changes the ball's trajectories, which re-rolls which of the
> demo's edge settings drive the loop into a near-singularity — the ball is
> purely rolling and the *mechanism* binds, with no hop involved at all.  That
> is [#22](https://github.com/caliburn-engineering/caliburn/issues/22)'s
> territory rather than this section's, and it is tracked as
> [#29](https://github.com/caliburn-engineering/caliburn/issues/29), which is a
> release blocker.
>
> It landed red on purpose.  The fix is independently sound — validated against
> a central difference on a smooth drive, where a difference *is* valid, and the
> two agree to four decimals — and it is the baseline every subsequent
> measurement has to be read against.  Holding it back is what entangled it with
> the restitution work in a single commit and left the demo's loss numbers
> unattributable to either.  The contract that Nominal holds every setting the
> sliders offer binds the ticket that *closes* #23, not every commit on the way
> there.
>
> One of the three is not a lost ball at all: Nominal now **saturates** the
> servos, which falsifies the `LqrPreset` blurb "no saturation" that
> `test_auto_balance` pins.  That is the blurb mechanism working as designed —
> a claim that stops being true fails a build rather than merely misleading
> somebody — and whether the answer is a re-measured blurb or a real regression
> is #29's to settle, not a string to quietly edit.

> **A separation is a claim about the plate's velocity, so it is only as good
> as that velocity.**
> The rates come through the velocity Jacobian, `-J_pose⁻¹ J_alpha`, which
> amplifies without bound near a singularity: an over-aggressive gain drives
> the plate to a Jacobian condition number of 2464, and the rates that come
> back would launch the ball a metre into the air off a 0.26 m/s nudge.  That
> hop is arithmetic, not physics.  `PlateMotion::rates_trustworthy` carries the
> credibility and the contact test declines to act without it, using the
> threshold the application already shows the user — condition 20, where its
> own readout turns red and says "Poor".  This is #22's lesson again: do not
> act on a solve you cannot validate.

> **An over-aggressive tuning does not merely overshoot — it throws the ball
> off the plate.**
> Q on ball position at 2000 loses the ball in 54 of 720 kick directions, with
> launches reaching 689 mm.  Before this the same tuning looked merely fast,
> because a glued ball cannot be thrown.  That is the honest ceiling on #19's
> aggressive preset.
>
> **These two figures were measured with the differenced estimator and have not
> been re-measured.**  Everything above about the shipped tuning has been; these
> have not, and by the same argument they are likely overstated.  Read them as
> the shape of the result rather than its size until
> [#29](https://github.com/caliburn-engineering/caliburn/issues/29) re-runs them.

### Trajectory tracking

A moving setpoint: circle, square, triangle, or the fixed point the loop has
always had.  Balancing at the centre proves stability; tracing a shape proves
*tracking*, and it is far more legible to someone who does not read a
pole-zero map.

The path writes `sp_x_mm_` / `sp_y_mm_` rather than going round them, so the
control law is untouched and the sliders become the readout of where the ball
is being sent — the same arrangement the servo sliders have under the loop.
Polygons are traversed at constant **speed**, not constant angle: a corner is
where the interesting behaviour is, and sweeping an angle would crawl through
it.

**Position on the path is a phase, and the phase is accumulated.**  `pathPoint`
and `pathVelocity` take a lap fraction in `[0, 1)`; `plate_view` holds it and
advances it by `dt / period_s` each step.  It used to be derived from the
simulation clock as `t / period_s`, which reads as a pure function of time and
is not one: changing the lap moves that quantity by `t dT / T^2`, and `t` is the
whole time the demo has been running.  Measured at 100 s, a nudge of the lap
slider from 10.0 to 9.5 s threw the setpoint **170 degrees** round the path, and
the loop hauled the ball across the plate after it.  The size slider did it too,
indirectly — a bigger path raises the lap floor, which pushes `period_s` up.

Accumulating instead means a lap change alters only the **rate**, from that
moment on, and a size change slides the setpoint **radially**: same angle,
bigger shape.  **A shape change re-seeds the phase** rather than carrying it,
because phase is not comparable across shapes: the circle's phase zero is at +x
and a polygon's first corner is at the top, so equal phase is a quarter of a lap
apart — measured, **169.7 mm** between a 120 mm circle and the square, at every
phase in the lap.  Carrying it threw the target to the far side of the path and
the loop hauled the ball across after it, into the workspace clip.  `phaseNearest`
picks the point on the new shape closest to where the setpoint already is, which
moves it by the least the two shapes allow: 35 mm circle to square, 60 mm circle
to triangle — exactly the `r/2` a triangle's inradius leaves — and nothing at all
when the shape does not change.

So all three trajectory controls now keep one promise: **the setpoint never
jumps.**  The lap changes its speed, the size slides it radially, and the shape
moves it as little as two different shapes permit.  Two rules go with it, and both live in `setpoint_path.h` rather
than in the panel so that the tests exercise them rather than their own copies:
`stepPath` owns the **order** — read the setpoint at the phase the frame opened
on, advance afterwards, so the position and the velocity handed to the loop are
the same instant and a polygon cannot straddle a corner between them — and
`clampPeriod` owns the **lap floor**, applied against the radius just dragged
rather than the one the last frame copied.  The reference *velocity* still steps on a lap change, which is
correct — the setpoint really was asked to travel faster — and
`test_setpoint_path` says so alongside the invariance.  All eight of that
file's original tests evaluated the path at fixed parameters, which is exactly
why none of them saw this; the ones added with the fix move a slider mid-run,
at four run lengths out to 1000 s, because the old error grew with `t` and a
brief test would have passed against it.

**Velocity feedforward, decided by measurement.**  The reference state has
always claimed the ball should be at the setpoint *and stationary*, which is
false the moment the setpoint moves — so the loop spent its effort fighting the
motion it was asked for.  Measured on a 120 mm circle at a ten-second lap:
**3.66 mm of mean error with it, 35.52 mm without** — nearly ten times.
Without it the ball does not follow the circle so much as sit inside it.

**It lives in the loop, as `BallReference`.**  It was applied for a while by
biasing the *measurement* handed to `legCommand`, from `plate_view` — identical
arithmetic, since only the difference enters `u = -K(x - x_ref)`, and cheaper
to write.  What made it wrong was everything around it: `auto_balance.h` went
on saying the regulator needs no feedforward while the application supplied one
two files away, and `test_trajectory` carried a third copy of it.  A second
implementation of the loop inside its own test harness is what let a real bug
hide once already, in #23.  So `legCommand` takes a reference *state* — where
the ball should be and how fast the setpoint is going — and there is one
implementation, in the file whose header describes the contract.

One thing deliberately left duplicated: **when** the reference velocity applies.
`plate_view` and both closed-loop harnesses each carry the same one-line rule —
zero it while the ball is airborne, because the plate cannot touch the ball and
`seen` is already the predicted landing point.  Collapsing it would mean either
`auto_balance` including `setpoint_path` or the reverse, and a loop that knows
about paths is a worse trade than a rule stated three times and documented once,
on `BallReference`.

**It has no gain of its own, and must not be given one.**  The reference
velocity is multiplied by K's velocity columns, which LQR designs from the
`x' ball` and `y' ball` weights — so the two sliders that own those numbers
already set how hard the loop chases a moving setpoint, and the LQR panel says
so beneath them.  A feedforward slider would be a second opinion about a number
K owns.  `test_auto_balance` pins both halves: that the reference velocity
enters as a velocity error, and that raising those two weights raises the
command a moving setpoint asks for.

The corner survives that: a square is tracked to about 14 mm at its corners
against a circle's 5.1 mm on the same size and lap — roughly three times.  That
is a bandwidth limit made visible, and it is the point of offering cornered
shapes at all.

**Three times is the claim; the square's own figure is not.**  It is a peak
sampled at 60 Hz exactly where the reference velocity steps, so which frame
lands nearest the corner decides it: a phase perturbation in the last bit of a
double moves it between 14.2 and 14.9 mm, while the circle's 5.072 mm holds to
a thousandth across the same change.  `test_trajectory` asserts the ratio for
that reason, and a fourth significant figure on the square would be a claim
about frame alignment rather than about the controller.

Decided in [#24](https://github.com/caliburn-engineering/caliburn/issues/24).

> **The sliders are bounded, and there are two bounds because there are two
> ways to ask for the impossible.**
> `kMaxSetpointSpeed` (250 mm/s) binds the large paths — beyond about 400 mm/s
> the plate loses the ball outright — and `kMinLapSeconds` (2 s) binds the
> small ones, where 250 mm/s round a 20 mm circle is a half-second lap and an
> angular rate the plate cannot follow however short the distance.  Neither
> implies the other.  A demo whose controls include a setting that breaks it is
> not offering a choice, it is offering a trap.

### Rolling-friction dead band

`RollingBallDynamics` opposes motion with `c_rr * (N/m) * sign(v)`, so a plate
tilted by less than `atan(c_rr) = 0.573 deg` **cannot start the ball moving at
all**: the tangential pull is `g sin(t)` and the resistance `c_rr g cos(t)`, so
the band closes where `tan(t) = c_rr`.

The derivation moved with #23 and the number did not.  Before the normal force
was computed, the resistance was scaled by `g` outright, which put the edge at
`asin(c_rr)` instead — 0.5729673 deg against 0.5729387 deg, a difference of
2.9e-5 deg, or one twenty-thousandth of the tolerance the test pins it to.  The
two agree to four decimals for any `c_rr` this small, so no test number changed;
what changed is which of them is *true*.  A state feedback has no integral term and therefore no way out: it
parks the ball wherever the tilt it is asking for falls inside the band.

This is not modelled in the linear plant at all, and it is why the LQR designer
does not open on unit weights.  Under `Q = R = I` the residual is 37 mm off
centre with the plate sitting at the dead-band edge — a loop that looks broken
while behaving exactly as designed.  The residual scales inversely with the
position gain, so the shipped default weights the ball position and lands
inside 2 mm.  Both the good case and the dead band itself are pinned by
`tests/test_auto_balance.cpp`; the number is not a memory.

> **Not "stiction" and not "the loop has steady-state error".**
> Stiction is a break-away force distinct from the sliding one, which this
> model does not have.  And a steady-state error suggests a gain that could be
> raised until it goes away — the band is a *region* of equilibria, and inside
> it the loop is not converging slowly, it is not moving.
