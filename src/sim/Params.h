#pragma once
#include <string>

namespace sim {

// Model parameters. Units: meters, seconds. Energies are measured in units where the rotational
// inertia chi = 1, so J, J4, T are in s^-2 and eta in s^-1. Speeds are handled internally as the
// dimensionless w = s / v0. See docs/SIMULATION_DESIGN.md section 3 for the equations.
struct Params {
    enum class Turning { Overdamped, Inertial };
    enum class SpeedControl { Linear, Marginal };
    enum class NeighborRule { Balanced, Nearest };

    int    N  = 300;
    double v0 = 12.0;               // reference cruise speed (m/s)

    // turning (heading) dynamics
    Turning turning = Turning::Inertial;
    double chi = 1.0;               // rotational inertia (sets the energy unit)
    double eta = 0.5;               // rotational friction
    double J   = 100.0;             // alignment strength
    double J4  = 0.0;               // quartic alignment strength (PRL 2026 model)
    double T   = 4.0;               // noise temperature

    // speed dynamics
    SpeedControl speedControl = SpeedControl::Linear;
    double g       = 700.0;         // linear speed-control stiffness, V = (g/2)(w-1)^2
    double lambda  = 0.7;           // marginal speed-control strength, V = lambda (w^2-1)^4
    double gammaS  = 140.0;         // speed friction
    double speedCoupling = 1.0;     // speed alignment = speedCoupling * J (1 = paper's unified model)

    // neighbors (topological)
    NeighborRule neighborRule = NeighborRule::Balanced;
    double balanceAngleDeg = 51.0;  // "balanced" rule: no two neighbors closer than this in angle
    int    nc = 7;                  // neighbor count for the plain nearest rule

    // cohesion forces (paper eq. G5, rescaled to meters) plus a hard core
    double rHardCore = 0.4;         // ~ starling wingspan (review 2018)
    double rEq       = 1.0;         // force-free distance ~ nearest-neighbor distance
    double rAttract  = 1.6;
    double cohesion  = 20.0;        // strength, same units as J
    double farAttract = 1.0;        // attraction beyond rAttract (1 in eq. G5)
    double coreStrength = 5.0;      // extra repulsion slope inside the hard core [ours]
    double cohesionSpeedGain = 1.0; // how strongly cohesion adjusts speed (fore-aft spacing) [ours]

    // integration
    double dt = 0.001;              // time step (s)
};

// The presets from docs/SIMULATION_DESIGN.md. Values are calibrated with the flocksim_cli tool;
// see docs for the calibration record.
enum class TurningPreset { Overdamped, Inertial, NonlinearInertial };
enum class SpeedPreset { Stiff, NearCritical, Marginal };

Params makePreset(TurningPreset turning, SpeedPreset speed, int N = 300);

const char *toString(TurningPreset p);
const char *toString(SpeedPreset p);

} // namespace sim
