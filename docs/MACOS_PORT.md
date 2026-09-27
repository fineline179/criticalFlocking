# macOS port notes (Cinder 0.9.2/Windows → Cinder 0.9.3/macOS)

This project was written in 2016–2017 in Visual Studio on Windows against Cinder 0.9.2. In September 2026
it was ported to build and run on macOS (Apple Silicon) against Cinder 0.9.3. This document records what
was changed, why, and what was left alone, so the context isn't lost.

Build instructions are in the [README](../README.md).

> **Later change (2026-09-27):** the simulation code this document describes (`include/swarming_spp/`,
> `Particle`, `ParticleController`) has since been replaced by a new simulation core in `src/sim/`. See
> [SIMULATION_DESIGN.md](SIMULATION_DESIGN.md). The macOS, Cinder, Retina and batched-drawing notes below
> still apply.

## Cinder version and setup

- **Cinder 0.9.3** (released Nov 2025) is the latest tagged release and is what the project targets.
  Cinder `master` (0.9.4dev) was deliberately *not* used: it is unreleased and replaces the macOS app
  layer with a GLFW-based rewrite.
- Cinder is expected in a `Cinder` directory next to this repo (`../Cinder`), overridable with
  `-DCINDER_PATH=...`. `scripts/setup_cinder.sh` clones it, patches it, and builds the static library
  with CMake (`Release` by default; `all` also builds `Debug`). The script is idempotent.
- **Cinder patch:** Cinder 0.9.3 does not compile with current macOS SDKs. Its bundled FreeType zlib
  (`src/freetype/gzip/ftzconf.h`) only typedefs `Byte` when `TARGET_OS_MAC` is undefined, expecting
  `MacTypes.h` to provide it, which recent SDKs no longer do. The script changes that guard to
  `#if !defined(Byte)`, which is the same fix Cinder `master` uses.
- Toolchain used: macOS 27, Apple clang 21 (Command Line Tools only, no Xcode), CMake 4.x from Homebrew.
  Cinder builds for `arm64;x86_64` and the app for the native arch. A libc++ warning about the 10.13
  deployment target is harmless.

## Build system

- A top-level `CMakeLists.txt` (using Cinder's `ci_make_app`) replaces the Visual Studio 2015
  `.sln`/`.vcxproj` files, which were deleted (they are in git history before this change). Those
  targeted Cinder 0.9.2 / Win32 / v140, and Cinder 0.9.3 dropped 32-bit Windows, so they could not
  have worked with the ported code anyway. The CMake project should in principle also work on Windows
  with Cinder 0.9.3, but that has not been tested.
- The build type defaults to `Release` (Cinder's CMake would default to `Debug`).
- `include/swarming_spp/hostile_environment.cpp` was never part of the VS project and is still not built.
- `.gitignore`: the old file ignored every `*.txt`, which would silently exclude `CMakeLists.txt`,
  so an exception was added. `build/`, `.DS_Store` and `imgui.ini` are now ignored. ImGui writes its
  window layout to `imgui.ini` in the current working directory.

## Code changes

### UI: AntTweakBar → Dear ImGui

Cinder 0.9.3 removed AntTweakBar (`cinder/params/Params.h`, `params::InterfaceGl`) and replaced it
with Dear ImGui (`cinder/CinderImGui.h`). The parameter panel was rebuilt in `FlockingApp::drawParamsWindow()`:

- The same parameters, ranges and read-only statistics (J, g, Temp, Balance Angle, Eye Distance,
  Bird From/To, Paused, Draw grid; Mean speed, Av Num Neighs, Polarization, Q_int).
- AntTweakBar's `keyIncr`/`keyDecr` bindings are reimplemented in `keyDown()` with the same keys and
  step sizes (see the README's controls table). `r` is new and resets the view.
- **Scene Rotation:** ImGui has no equivalent of AntTweakBar's quaternion rotation widget. It is
  replaced by left-drag in the scene to orbit the camera around the flock center, mouse wheel to zoom
  (eye distance), and a "Reset view" button. Mouse events over the ImGui panel are consumed by ImGui
  and don't rotate the scene.
- `birdFrom` / `birdTo` changed from `float` to `int`, and their slider maximum is now `mN - 1`
  instead of a hardcoded 9, which resolves the old TODO.
- `mQ_int` is now initialized to 0 (it used to display garbage until the simulation was unpaused).
- `resize()` now updates the camera aspect ratio.

### Performance: batched drawing (~6 fps → 60 fps)

On macOS the app initially ran at ~6 fps. Profiling showed the simulation itself was negligible; the
time went into Cinder's immediate-mode draw helpers (`gl::drawCube`, `gl::drawLine`, `gl::drawSolidRect`,
`gl::drawVector`, …). Each one uploads its vertices into a single shared VBO with `glBufferSubData`,
and on Apple's Metal-backed OpenGL that forces a CPU/GPU sync on every call. With 300 birds plus ~730
grid lines per frame, that capped the frame rate.

- Birds are drawn with a cached unit-sphere `gl::Batch` plus a model matrix (`Particle::draw` now takes
  the batch; `FlockingApp::mSphereBatch`). This is equivalent to `gl::drawSphere(mPos, radScale * mRadius, 8)`.
  See "Bird highlighting" below; the first version of the port used a cube batch here.
  `ParticleController` (unused by the app) was updated to match.
- The grid is drawn as a single `gl::begin(GL_LINES)` … `gl::end()` batch instead of one `gl::drawLine`
  per segment.
- The C_sp(r) bar graph is drawn as one `GL_TRIANGLES` batch, and its axes as one `GL_LINES` batch,
  instead of per-bar `gl::drawSolidRect` and per-axis `gl::drawLine` calls.
- The neighbor arrows of a highlighted bird are drawn as one `GL_LINES` batch, with arrowheads from a
  cached unit-cone batch (`FlockingApp::drawArrowHead`, same geometry as `gl::drawVector`'s head).
  With these still per-call, turning on highlighting dropped the frame rate to 30 fps.
- Only the bounding cube and the flock-velocity arrow still use per-call helpers (3 calls per frame).
  If more drawing is added, batch it; the cost per call grows with how much the GPU has queued, so a few
  extra calls can push a frame over budget.

### Retina (high-density) display

On this macOS version the `NSOpenGLView` gets a Retina (2×) drawable even when Cinder's high-density
mode is off. Cinder's viewport matched that drawable, but Cinder reported a content scale of 1, so ImGui
rendered at half size in the lower-left corner with mismatched mouse coordinates. The fix:

- `settings->setHighDensityDisplayEnabled(true)` in `prepareSettings`, so Cinder, the drawable and
  ImGui all agree on 2×. Do **not** try to force a 1× viewport; that is the wrong fix.
- ImGui works in framebuffer pixels, so its style and font are scaled by `getWindowContentScale()`
  (`ScaleAllSizes`, `FontScaleDpi`), and the panel auto-sizes to its contents.

### Minor

- `M_PI` definitions in `swarming_spp/behavior.cpp` and `community.cpp` are guarded with `#ifndef M_PI`,
  since macOS's `<math.h>` already defines it.

### Bird highlighting (fixed after the initial port)

In the original code, birds were drawn as cubes because `gl::drawSphere` drew stray lines from the
origin in Cinder 0.9.2 (the old TODO in `Particle::draw`). Cubes ignored `radScale`, so the highlighted
"from" (white) and "to" (blue) birds weren't drawn at their intended 1.6× size. Birds are now spheres of
radius `radScale * mRadius` (0.05, so slightly larger than the old 0.07-wide cubes). The highlighted birds
are therefore 1.6× bigger, and the stray-lines problem doesn't occur with the batch approach.

## Known issues left alone (pre-existing)

- `Particle::pullToCenter` uses `dirToCenter.length()`, which for a glm vector returns the component
  count (3), not the magnitude. It should be `glm::length(dirToCenter)`. The app never calls it.
- `gl::rotate(mSceneRotation)` at the top of `FlockingApp::update()` has no effect, because the
  matrices are replaced by `gl::setMatrices(mCam)` right after.

## How the port was verified

The terminal session didn't have macOS Screen Recording permission, so `screencapture` couldn't be used.
Instead, a throwaway copy of the app (not committed) hooked the window's post-draw signal. It saved
frames with `copyWindowSurface()`, injected synthetic input through `getWindow()->emitMouseDown/Drag/Wheel`
and `emitKeyDown` (so ImGui gets first claim on events, as with real input), logged state, and quit.
(Screen Recording later turned out to work from iTerm2. Capturing just the app window, even when it is
covered, works with `screencapture -o -l <CGWindowID>`, where the window ID comes from
`CGWindowListCopyWindowInfo` using the app's entry with nonzero bounds.)
That confirmed:

- The flock renders centered on its center of mass, at 60 fps, in both the paused and running states.
- The flock statistics look right (polarization ≈ 0.99, ~7 neighbors for n_c = 8, mean speed ≈ v0).
- Bird highlighting shows the neighbor arrows.
- Dragging over the scene rotates it and dragging over the ImGui panel doesn't; the wheel zooms.
- Every keyboard shortcut changes the right parameter by the right step.
- Both `Release` and `Debug` configurations build, and `scripts/setup_cinder.sh` works from a fresh
  clone and on an existing checkout.
