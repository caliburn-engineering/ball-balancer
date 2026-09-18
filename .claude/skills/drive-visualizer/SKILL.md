---
name: drive-visualizer
description: Launch, focus, click and screenshot the Ball-Balancer visualizer GUI on Hyprland/Wayland. Use when verifying a UI change in the app, taking screenshots of panels, or driving the app to a specific state.
---

# Driving the Ball-Balancer visualizer

Desktop GLFW/OpenGL app. On this machine (Hyprland + Wayland) it is driven with
`hyprctl` for focus, `ydotool` for input, and `grim` for screenshots.

## Launch — switch to workspace 3 FIRST

**Workspace 3 is the human's designated Caliburn screen.** Always switch there
before launching or driving anything, and keep the mouse there. The switch is
the human's signal that the agent has taken the pointer and they should stop
interfering; driving the app on whatever workspace happens to be visible means
fighting them for the mouse and screenshotting their other work.

```bash
hyprctl dispatch "hl.dsp.send_shortcut({mods='SUPER', key='3'})"   # literally super+3
sleep 0.5
hyprctl activeworkspace -j | python3 -c "import json,sys; print(json.load(sys.stdin)['id'])"
```

Then launch from the repo root — the app reads `imgui.ini` and `vendor/fonts/`
relative to cwd, and opens on the *active* workspace, so the order matters.

```bash
cd projects/ball-balancer
setsid ./build/visualizer >/tmp/viz.log 2>&1 &
sleep 4
```

If the app is already running somewhere else, pull it over rather than
switching away:

```bash
hyprctl dispatch "hl.dsp.focus({window='address:$A'})"
hyprctl dispatch "hl.dsp.window.move({workspace='3'})"
```

The window title is `Ball-Balancer — Caliburn` since the merge (issue #13);
`Linear System Analyzer` no longer matches anything.

Finally park the pointer on workspace 3 (see **Click**) so the human can see at
a glance that the agent is driving.

## Find the window

```bash
A=$(hyprctl clients -j | python3 -c "
import json,sys
for c in json.load(sys.stdin):
    if 'Ball-Balancer' in c.get('title',''):
        print(c['address'])")
```

Geometry (`at` and `size`) comes from the same JSON — you need it to convert
screenshot pixel coordinates to screen coordinates.

## Focus it — REQUIRED before any click

**This Hyprland build uses a Lua dispatch API.** The old string form
(`hyprctl dispatch focuswindow address:0x...`) fails with a Lua parse error, and
so does `hl.dsp.window.focus` — there is no `focus` under `hl.dsp.window`.
The dispatcher is `hl.dsp.focus`, and it takes a **table**:

```bash
hyprctl dispatch "hl.dsp.focus({window='address:$A'})"
```

Without this, `ydotool` clicks go to whatever window currently has focus —
usually the terminal — and the app silently ignores them.

Discover any other dispatcher the same way, rather than guessing:

```bash
hyprctl repl 'local t={} for k,v in pairs(hl.dsp) do t[#t+1]=k end table.sort(t) return table.concat(t,", ")'
hyprctl dispatch "hl.dsp.focus({})"   # error message lists the accepted keys
```

`hyprctl repl <lua>` evaluates and prints; `hyprctl dispatch <lua>` wraps the
argument in `hl.dispatch(...)`.

## Click

`ydotoold` must be running (`pgrep ydotoold`). Coordinates are **screen**
coordinates: window `at` + pixel offset within the screenshot.

```bash
hyprctl dispatch "hl.dsp.cursor.move({x=<X>, y=<Y>})"; sleep 0.4
hyprctl dispatch "hl.dsp.cursor.move({x=<X>, y=<Y>})"; sleep 0.4
hyprctl cursorpos          # confirm the warp landed before committing to a click
ydotool click 0xC0
```

`0xC0` is press+release of the left button. Re-focus the window before each
click burst if you ran anything that may have stolen focus.

**Warp twice.** A single warp followed immediately by a click has the press
arrive while ImGui still believes the pointer is where it was, and the release
at the new position turns the whole thing into a *drag* — which silently
rearranges the dock layout, moves floating windows, and edits whatever slider
was under the old position. The symptom is a click that "did nothing" plus a
layout that quietly changed. The second warp gives ImGui a frame to catch up.

`ydotool mousemove -x <dx> -y <dy>` (relative) is no substitute: pointer
acceleration scales the delta, so the cursor does not land where the arithmetic
says it will.

## Screenshot

```bash
grim -g "<x>,<y> <w>x<h>" /tmp/shot.png     # region
```

Crop tightly to the panel you care about — a full 1896x1030 frame is ~1.9 MB and
mostly plots you are not checking. Then read the PNG.

## Gotchas

- **A blank or unchanged screenshot after a click usually means focus, not a
  code bug.** Re-check `hyprctl activewindow -j` before concluding anything.
- The window may be tiled small. It resizes when it is the only tiled client.
- **A hollow box has two different causes, and they need different fixes.**
  Measured against the bundled `NotoSans-Regular.ttf` (issue #15):

  | Block | In the font | Verdict |
  |---|---|---|
  | Arrows U+2190-21FF | 0 / 112 | Write ASCII: `<-`, `->` |
  | Mathematical Operators U+2200-22FF | 1 / 256 (only U+2212 MINUS) | Write ASCII: `inf`, `arg`, `~=` |
  | General Punctuation U+2000-206F | 111 / 112 | Fine — the range is requested |
  | Phonetic Extensions U+1D00-1D7F | 128 / 128 | Fine — the range is requested |
  | Greek, super/subscripts, Latin-1 | covered | Fine |

  So: a glyph the **font** lacks can never be made to draw, and must be
  rewritten in ASCII. A glyph the font has but `glyph_ranges` in
  `visualizer.cpp` does not **request** boxes just as convincingly, and is
  fixed by adding the range. Em dash and bullet were the second kind and drew
  as boxes in 15 UI strings until #15.

  To check a candidate character before using it, read the TTF's `cmap`
  rather than guessing.
- Watch for **C++ hex-escape greediness** in UI strings: `"\xcf\x84f"` is parsed
  as `\xcf` then `\x84f` (out of range), not `tau` + `f`. Split the literal:
  `"\xcf\x84" "f"`.

## Click discipline

**Screenshot the region first, locate the widget, then click.** The model panel
scrolls, so a coordinate that was a combo box a moment ago can be a matrix
slider now — and a stray click on a `SliderFloat` silently edits the plant.
If a value looks wrong afterwards, restart the app rather than trying to undo.

Scrolling: `ydotool mousemove -w -x 0 -y -12` scrolls **down**, positive `-y`
scrolls up. `ydotool mousemove -a` does **not** work on this machine — the
cursor does not move — which is why focus uses `hl.dsp.cursor.move`.

Combos near the bottom of the window flip their popup upward; move the panel so
the widget sits mid-height before clicking, or the entry you want will not be
where you predict.

**Never use `pkill -f <pattern>` where the pattern appears in your own command
line** — `-f` matches the agent's shell too, and the session dies mid-task.
This has happened twice: `pkill -f build/visualizer` and
`pkill -f "http.server 8731"`. Use an exact process name, or resolve the PID:

```bash
pkill -x visualizer                                    # exact name
PID=$(ss -lptn 'sport = :8731' | grep -oP 'pid=\K[0-9]+' | head -1)   # by port
```

## Workspaces and scroll traps

- **The window may be on another workspace.** `grim` captures the *visible*
  output and `ydotool` clicks land on the *visible* workspace, so a stale
  window elsewhere means you screenshot someone else's app and click into it —
  this has already put stray clicks into the human's browser once. Switch to
  workspace 3 first, confirm with `hyprctl activeworkspace -j`, and confirm the
  focused title with `hyprctl activewindow -j` before the first click.
  `hl.dsp.focus({window=...})` also switches workspace as a side effect.
- **ImPlot eats the scroll wheel to zoom.** Scrolling with the cursor over any
  plot zooms that plot instead of scrolling its panel, and there is no undo —
  restart the app. Scroll panels from a non-plot column (the model panel is
  safe at x ~170).
