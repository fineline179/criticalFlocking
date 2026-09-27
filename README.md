# _criticalFlocking

A real-time starling-flock simulation built on the statistical model of
[Bialek et al. 2013 (arXiv:1307.5563)](https://arxiv.org/abs/1307.5563), extended with inertial turning
dynamics, alternative speed controls and a predator. Its purpose is to show when a threat sensed by a few
birds spreads through the whole flock "all as one", and when it stays local.

The original version (a reproduction of the paper's eqs. G1–G5, based on
[Swarming SPP](https://github.com/david-mateo/swarming-spp)) is on the `master` branch.

## Documentation

- [docs/CRITICALITY_AND_DYNAMICS.md](docs/CRITICALITY_AND_DYNAMICS.md): what the paper establishes (static
  statistics) vs. what a real-time predator response needs (the dynamics).
- [docs/SIMULATION_DESIGN.md](docs/SIMULATION_DESIGN.md): the model, its presets, the predator, how it was
  calibrated against published starling data, and what the results do and don't show.
- [docs/MACOS_PORT.md](docs/MACOS_PORT.md): the port from Visual Studio / Cinder 0.9.2 to macOS / Cinder 0.9.3.

## Building (macOS)

Requires the Xcode command line tools (`xcode-select --install`), CMake (`brew install cmake`) and
[Cinder](https://libcinder.org) 0.9.3.

1. Fetch and build Cinder. This clones Cinder 0.9.3 into `../Cinder` (next to this repo), applies a
   one-line fix needed for recent macOS SDKs, and builds it (pass `all` to also build the Debug library):

   ```sh
   scripts/setup_cinder.sh
   ```

   To use a Cinder checkout somewhere else, set `CINDER_PATH=/path/to/Cinder` for the script and pass
   `-DCINDER_PATH=/path/to/Cinder` to CMake below.

2. Build and run the app:

   ```sh
   cmake -S . -B build
   cmake --build build -j
   open build/Release/criticalFlocking/criticalFlocking.app
   ```

   Rebuild (`cmake --build build -j`) after any code change before opening the app. For a debug build
   use `cmake -S . -B build/debug -DCMAKE_BUILD_TYPE=Debug`.

## Using the app

- **Flock panel:** choose the turning dynamics (A overdamped, B inertial, C nonlinear inertial) and the
  speed control (S1 stiff, S2 near-critical, S3 marginal). Hover a preset for its source and status. Also:
  time scale (slow motion), flock size, color mode, and a bird whose neighbor links to show.
- **Predator panel:** launch an attack, or turn on automatic attacks; set the detection radius, how hard
  detecting birds turn away and speed up, and the rate of false alarms. After each attack it reports
  what happened. These are model predictions, not facts about real starlings.
- **Validation panel:** live comparison of the simulated flock with published measurements of real
  starling flocks (pass/fail per row; hover a row for its source).
- **Graphs (lower left):** heading (green) and speed (blue) correlation vs. distance; the yellow line
  marks the correlation length.
- **Colors:** by default, each bird's color shows how far it has turned since the last attack began
  (white → red at 60°), so a spreading turn shows up as a spreading red front.

| key     | action                              |
|---------|-------------------------------------|
| space   | pause / run                         |
| a       | predator attack                     |
| f       | toggle false alarms                 |
| 1 2 3   | turning preset A / B / C            |
| 4 5 6   | speed preset S1 / S2 / S3           |
| [ ]     | slower / faster time                |
| c       | cycle color mode                    |
| n       | new flock                           |
| g       | toggle grid                         |
| r       | reset view                          |
| drag / wheel | rotate / zoom                  |

## Headless tool

`build/flocksim_cli` runs the same simulation without graphics, for calibration and validation:

```sh
build/flocksim_cli stats  --turning inertial --speed marginal       # snapshot statistics
build/flocksim_cli turn   --turning overdamped --pulse 0.3          # does a brief turn spread?
build/flocksim_cli attack --turning inertial --seed 3               # one predator attack
build/flocksim_cli speedpush --speed stiff                          # does a speed-up spread?
build/flocksim_cli probe  --turning nonlinear                       # spontaneous fluctuations
```

Any parameter can be overridden with `--set name=value` (see `tools/flocksim_cli.cpp`).

## Code layout

- `src/sim/`: the simulation core (no Cinder dependency): `Flock` (state, neighbors, integrators),
  `Params` (parameters and presets), `Predator`, `Measure` (statistics, turn-front tracking,
  fluctuation probe).
- `src/FlockingApp.cpp`: the Cinder app (rendering, UI, validation panel).
- `tools/flocksim_cli.cpp`: the headless driver.
