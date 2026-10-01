# Ball-Balancer

A 3-RRS table balancing a ball under real rolling dynamics, and the analyzer
that models a plant, closes a loop around it, and shows the consequences across
Bode, Nyquist, pole-zero and step-response views. One executable over two
halves.

Named `linear-analyzer` until
[#20](https://github.com/caliburn-engineering/caliburn/issues/20): it is the
survivor of the merge and carries the whole history. The pre-merge
ball-balancer repository stands unmodified beside it as
`projects/ball-balancer-legacy`, as a historical record.

## CONTEXT.md is the vocabulary

`CONTEXT.md` holds this project's ubiquitous language: **loop**, **RGA**,
**Channel Share**, **structurally dead channel**, **pairing grid**, **plate
view**, **balance loop**, **leg command vs leg angle**, **assembly mode**,
**workspace vs servo box**, **contact**, **fillet**. Read the entry before
using a term in code, in a comment, or anywhere on screen — several entries
record a plausible synonym that is *wrong* here together with the reason, and
drifting to it reintroduces a resolved confusion.

Decisions live inline in that glossary with their issue link. Record a new one
the same way when it attaches to a term, and in `docs/adr/` when it does not.
`docs/plans`, `docs/specs` and `docs/research` are a dated archive from an
earlier workflow; read them for history and write new decisions in the two
places above.

## One transcription

A fact about the plant, or the order of operations inside a frame, lives in
exactly one place. Five copies of the frame loop cost three bugs before
[#30](https://github.com/caliburn-engineering/caliburn/issues/30): the
harnesses read the servo rate before the step where the application read it
after — a factor of 1.40 at 60 Hz against a 0.05 s lag — so a green suite was
measuring a plate that was not the shipped one.

| The one place | What it owns |
|---|---|
| `stepSim` (`src/sim_step.h`) | What a frame IS: setpoint, loop, servos, pose, plate motion, ball, in that order |
| `cascadePlate` / `cascadeDesign` | The shipped plant, assembled from a parameter list |
| `SimPlate::feasible` | The corner fillet, sized by the mechanism |
| `kPlateBall` (`src/ball_sim.h`) | The simulated ball |
| `tests/cascade_fixture.h` | The preset, gain and ball every closed-loop test measures |

A harness decides what to measure and when to shove the ball; `stepSim` decides
what a frame is. `test_ball_contact` is the deliberate exception: it drives
`plateMotion` and `stepBallContact` against plates built by hand, which is the
layer `stepSim` is built ON.

## Both targets

Every change builds for desktop (GLFW + glad + OpenGL 3.3 core) and for the web
(Emscripten + WebGL2). The renderer draws through the GL subset shared by both,
so a shader carries one `#version` line per target ahead of one shared body,
and `glLineWidth` above 1.0 is decorative under WebGL2 — every element of the
scene is told apart by colour.

Emscripten supplies GLFW and the GL headers through its own ports, so a new
desktop dependency breaks the web link silently. Build both before calling a
change done.

```sh
# desktop
cmake -S . -B build && cmake --build build -j && ctest --test-dir build

# web
. ~/emsdk/emsdk_env.sh
emcmake cmake -S . -B build-web && cmake --build build-web -j
```

Tests are desktop-only. Both build directories are gitignored.

## The frame contract

`PlateView` owns no GLFW window, no ImGui context and no main loop —
`visualizer.cpp` owns all three, and reconciling that was the whole of the
merge. The caller fixes the order: `attach()` before the ImGui backend is
initialised, `initGL()` once the GL loader is up, then `step()` →
`drawPanels()` → `drawScene()` every frame.

`handDesignToPlate` hands the plate the analyzer's current design once per
frame, and that is the only traffic between the halves: the plate never writes
`AppState`.

The 3D scene shows through the **central dock node**, which stays empty. A
window docked there paints over the plate.

## Before writing a numerical routine

Look in `../../reference/` and `../../knowledge/` first, and say what you
found. `reference/` is the workspace's golden source, and it is downstream of
the first project that needs a thing and upstream of every one after. This
project carries its own LQR, PID and servo lag because nobody looked;
reconciling them with `reference/` is open work.

Golden source is copied, never built against (caliburn ADR-0008). Take what
you need — a whole file into `src/`, or an excerpt into the source file it
fits — and put an origin comment on it:
`// From caliburn reference/<path> @ <short commit>`. No build file here may
point at `reference/` or anything else outside this repository: a lone clone
must configure, build and pass `ctest`.
