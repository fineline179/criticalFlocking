# _criticalFlocking

This is a reproduction of eqs. G1-G5 in [Bialek - 1307.5563](https://arxiv.org/abs/1307.5563)

requires [Cinder](https://libcinder.org) 0.9.3

Based on [Swarming SPP](https://github.com/david-mateo/swarming-spp), which has been (substantially) 
 modified in order to implement the balanced nearest-neighbor algorithm described in 1307.5563.

## Building (macOS)

The project was originally built with Visual Studio on Windows against Cinder 0.9.2. See
[docs/MACOS_PORT.md](docs/MACOS_PORT.md) for what changed in the port to macOS / Cinder 0.9.3 and why.

Requires the Xcode command line tools (`xcode-select --install`) and CMake (`brew install cmake`).

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

   For a debug build use `cmake -S . -B build/debug -DCMAKE_BUILD_TYPE=Debug`.

## Controls

The simulation starts paused; press space (or untick "Paused") to run it.

- Drag in the scene to rotate the view, scroll to zoom.
- The "Flocking" panel sets the model parameters and shows flock statistics. The bar graph in the lower
  left is the speed correlation C_sp(r) as a function of bird separation.
- Keyboard shortcuts:

  | key     | action                    |
  |---------|---------------------------|
  | space   | pause / run               |
  | f / d   | J +/-                     |
  | v / c   | g +/-                     |
  | y / t   | Temp +/-                  |
  | o / i   | Balance Angle +/-         |
  | s / w   | Eye Distance +/-          |
  | . / ,   | Bird From +/- (-1 = none) |
  | l / k   | Bird To +/-               |
  | g       | toggle grid               |
  | r       | reset view                |
