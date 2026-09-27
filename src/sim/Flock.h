#pragma once
#include "Params.h"
#include "Vec3.h"
#include <cstdint>
#include <random>
#include <utility>
#include <vector>

namespace sim {

// A flock of N birds. Each bird has a position x, unit heading u, spin sigma (perpendicular to u;
// turn rate |sigma|/chi) and dimensionless speed w = s/v0.
//
// Heading dynamics (docs/SIMULATION_DESIGN.md 3.3):
//   Inertial:    du/dt = (1/chi) sigma x u,   dsigma/dt = u x F - (eta/chi) sigma + u x xi
//   Overdamped:  eta du/dt = P_perp(u) F + noise
// Speed dynamics (3.4):
//   gammaS dw/dt = -dH/dw + 2T/(N w_mean) + (cohesion drive) + b + noise
// where F = -dH/du plus cohesion and the external heading field h (e.g. from a predator).
class Flock {
  public:
    Flock(const Params &params, uint64_t seed);

    // Place birds in a ball at roughly the nearest-neighbor spacing, all heading along +x.
    void reset();

    // Advance by `seconds` of simulated time in steps of params.dt, recomputing neighbors every
    // `neighborInterval` seconds.
    void advance(double seconds, double neighborInterval = 1.0 / 60.0);
    void step();
    void updateNeighbors();

    const Params &params() const { return mParams; }
    // Parameters that are safe to change on a running flock (not N)
    void setParams(const Params &params);

    int    size() const { return (int) x.size(); }
    double time() const { return mTime; }
    double speed(int i) const { return w[i] * mParams.v0; }
    Vec3   velocity(int i) const { return u[i] * (w[i] * mParams.v0); }
    // turn rate (rad/s) from the last step
    double turnRate(int i) const { return mTurnRate[i]; }

    // State (public for measurement and rendering)
    std::vector<Vec3>   x, u, sigma;
    std::vector<double> w;

    // External fields, reset to zero by the caller (e.g. the predator) as needed
    std::vector<Vec3>   h;   // heading field
    std::vector<double> b;   // speed field

    // Directed topological neighbor lists (who bird i listens to) and the symmetrized
    // interaction weights n_ij = (nhat_ij + nhat_ji)/2 used in the energy.
    std::vector<std::vector<int>>                   directed;
    std::vector<std::vector<std::pair<int, double>>> coupled;

  private:
    void computeForces(std::vector<Vec3> &F, std::vector<double> &Fw);
    void stepInertial();
    void stepOverdamped();
    void stepSpeedAndPosition(const std::vector<double> &Fw);
    double gauss() { return mNormal(mRng); }
    Vec3   gaussVec() { return { gauss(), gauss(), gauss() }; }

    Params                           mParams;
    std::mt19937_64                  mRng;
    std::normal_distribution<double> mNormal{ 0.0, 1.0 };
    double                           mTime = 0;
    double                           mSinceNeighbors = 1e9;
    std::vector<Vec3>                mF, mUprev;
    std::vector<double>              mFw, mTurnRate;
};

} // namespace sim
