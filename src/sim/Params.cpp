#include "Params.h"

namespace sim {

// Preset values, calibrated with tools/flocksim_cli against the targets in
// docs/SIMULATION_DESIGN.md section 5 (calibration record in section 10 there).
Params makePreset(TurningPreset turning, SpeedPreset speed, int N)
{
    Params p;
    p.N = N;

    // Shared by all presets, so that snapshot statistics match (same energy, same noise)
    p.v0                = 12.0;    // m/s (review 2018: 12 +- 2.5)
    p.chi               = 1.0;
    p.J                 = 200.0;   // with B: turn fronts ~14-24 m/s (data: 20-40 m/s)
    p.T                 = 20.0;    // polarization ~0.97-0.98 (data: 0.96 +- 0.03)
    p.gammaS            = 280.0;   // speed relaxation J n_c / gammaS ~ 5 /s [ours]
    p.speedCoupling     = 1.0;     // unified model (Bialek 2014)
    p.rHardCore         = 0.4;
    p.rEq               = 1.3;     // with farAttract/coreStrength: r1 ~0.7-1.0 m (data 0.68-1.51)
    p.rAttract          = 1.9;
    p.farAttract        = 0.15;
    p.coreStrength      = 10.0;    // stiffer cores heat the inertial flock (steering instability)
    p.cohesionSpeedGain = 30.0;    // fore-aft spacing held by speed changes [ours]

    switch (turning) {
        case TurningPreset::Overdamped:
            p.turning  = Params::Turning::Overdamped;
            p.eta      = 37.0;     // sqrt(J n_c chi): same local response time as B
            p.cohesion = 320.0;    // calibrated so the resting flock matches B's size and spacing
            break;
        case TurningPreset::Inertial:
            p.turning  = Params::Turning::Inertial;
            p.eta      = 2.0;      // damping time 2 chi/eta = 1 s > crossing time; >= 1 needed
                                   // to damp the cohesion/steering instability
            p.cohesion = 40.0;
            break;
        case TurningPreset::NonlinearInertial:
            p.turning  = Params::Turning::Inertial;
            p.eta      = 60.0;     // heavily damped: spontaneous fluctuations overdamped
            p.J4       = 30000.0;  // quartic stiffening for large misalignments (PRL 2026)
            p.cohesion = 110.0;
            p.dt       = 0.0005;   // quartic term is stiff
            break;
    }

    switch (speed) {
        case SpeedPreset::Stiff:
            p.speedControl = Params::SpeedControl::Linear;
            p.g            = 1.0 * p.J * p.nc;     // g / (J n_c) = 1: speed changes stay local
            break;
        case SpeedPreset::NearCritical:
            p.speedControl = Params::SpeedControl::Linear;
            p.g            = 1e-3 * p.J * p.nc;    // g / (J n_c) = 1e-3 (Bialek 2014)
            break;
        case SpeedPreset::Marginal:
            p.speedControl = Params::SpeedControl::Marginal;
            p.lambda       = 7e-3 * p.J;           // lambda / J as in Cavagna 2022's simulations
            break;
    }
    return p;
}

const char *toString(TurningPreset p)
{
    switch (p) {
        case TurningPreset::Overdamped: return "A. Overdamped";
        case TurningPreset::Inertial: return "B. Inertial";
        case TurningPreset::NonlinearInertial: return "C. Nonlinear inertial (experimental)";
    }
    return "";
}

const char *toString(SpeedPreset p)
{
    switch (p) {
        case SpeedPreset::Stiff: return "S1. Stiff";
        case SpeedPreset::NearCritical: return "S2. Near-critical (linear)";
        case SpeedPreset::Marginal: return "S3. Marginal";
    }
    return "";
}

} // namespace sim
