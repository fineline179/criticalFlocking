#include "Flock.h"
#include <algorithm>
#include <cmath>

namespace sim {

namespace {
const double kPi = 3.14159265358979323846;
// candidates examined (nearest first) when building the balanced neighbor list
const int kBalancedCandidates = 60;
const int kMaxNeighbors       = 16;
} // namespace

Flock::Flock(const Params &params, uint64_t seed) : mParams(params), mRng(seed)
{
    reset();
}

void Flock::setParams(const Params &params)
{
    int N     = mParams.N;
    mParams   = params;
    mParams.N = N;
    if (mParams.turning == Params::Turning::Overdamped)
        for (auto &s : sigma) s = Vec3();
    mSinceNeighbors = 1e9;
}

void Flock::reset()
{
    const int N = mParams.N;
    x.assign(N, Vec3());
    u.assign(N, Vec3(1, 0, 0));
    sigma.assign(N, Vec3());
    w.assign(N, 1.0);
    h.assign(N, Vec3());
    b.assign(N, 0.0);
    mF.assign(N, Vec3());
    mFw.assign(N, 0.0);
    mUprev = u;
    mTurnRate.assign(N, 0.0);
    mTime = 0;

    // random positions in a ball, at density ~0.2 / r_eq^3 (review 2018: rho r1^3 = 0.14-0.25),
    // rejecting points closer than the hard core
    double                                 radius = mParams.rEq * std::cbrt(3.0 * N / (4.0 * kPi * 0.2));
    std::uniform_real_distribution<double> uni(-1.0, 1.0);
    for (int i = 0; i < N; i++) {
        for (int attempt = 0; attempt < 1000; attempt++) {
            Vec3 p(uni(mRng), uni(mRng), uni(mRng));
            if (length2(p) > 1.0) continue;
            p *= radius;
            bool ok = true;
            for (int j = 0; j < i && ok; j++)
                ok = length2(x[j] - p) > mParams.rHardCore * mParams.rHardCore * 2.25;
            if (ok) { x[i] = p; break; }
        }
    }
    mSinceNeighbors = 1e9;
    updateNeighbors();
}

void Flock::updateNeighbors()
{
    const int N = size();
    directed.assign(N, {});
    std::vector<std::pair<double, int>> cand(N);

    const double cosMu = std::cos(mParams.balanceAngleDeg * kPi / 180.0);
    for (int i = 0; i < N; i++) {
        int m = 0;
        for (int j = 0; j < N; j++)
            if (j != i) cand[m++] = { length2(x[j] - x[i]), j };

        if (mParams.neighborRule == Params::NeighborRule::Nearest) {
            int k = std::min(mParams.nc, m);
            std::partial_sort(cand.begin(), cand.begin() + k, cand.begin() + m);
            for (int a = 0; a < k; a++) directed[i].push_back(cand[a].second);
        }
        else {
            // Balanced: take candidates nearest-first, rejecting any that lie within the balance
            // angle of an already-accepted (closer) neighbor.
            int k = std::min(kBalancedCandidates, m);
            std::partial_sort(cand.begin(), cand.begin() + k, cand.begin() + m);
            std::vector<Vec3> dirs;
            for (int a = 0; a < k && (int) dirs.size() < kMaxNeighbors; a++) {
                Vec3 d        = normalize(x[cand[a].second] - x[i]);
                bool shadowed = false;
                for (const Vec3 &e : dirs)
                    if (dot(d, e) > cosMu) { shadowed = true; break; }
                if (!shadowed) {
                    dirs.push_back(d);
                    directed[i].push_back(cand[a].second);
                }
            }
        }
    }

    // symmetrized weights n_ij = (nhat_ij + nhat_ji) / 2
    coupled.assign(N, {});
    auto add = [this](int i, int j) {
        for (auto &e : coupled[i])
            if (e.first == j) { e.second += 0.5; return; }
        coupled[i].push_back({ j, 0.5 });
    };
    for (int i = 0; i < N; i++)
        for (int j : directed[i]) {
            add(i, j);
            add(j, i);
        }
    mSinceNeighbors = 0;
}

void Flock::computeForces(std::vector<Vec3> &F, std::vector<double> &Fw)
{
    const Params &p = mParams;
    const int     N = size();
    const double  Js = p.J * p.speedCoupling;
    const double  coreForce = 0.25 * (p.rHardCore - p.rEq) / (p.rAttract - p.rHardCore);

    // Volume ("entropic") factor for speed. In the full-velocity models this is s^2 for the
    // flock's mean velocity only (each bird's own factor cancels against its sideways velocity
    // fluctuations), giving a push 2T / w_mean on the group speed, shared among the N birds
    // (Cavagna et al. 2022). It is what inflates group speed under weak speed control.
    double wMean = 0;
    for (double wi : w) wMean += wi;
    wMean /= N;
    const double entropic = 2.0 * p.T / (N * wMean);

    for (int i = 0; i < N; i++) {
        // alignment (quadratic + quartic) and speed alignment, from the energy H
        Vec3   Fi;
        double Fwi = 0;
        for (const auto &[j, n] : coupled[i]) {
            double c = dot(u[i], u[j]);
            Fi += u[j] * (n * (p.J + 4.0 * p.J4 * (1.0 - c)));
            Fwi -= n * Js * (w[i] - w[j]);
        }

        // cohesion (paper eq. G5 with a hard core); not derived from H
        Vec3 coh;
        if (!directed[i].empty()) {
            for (int j : directed[i]) {
                Vec3   d   = x[j] - x[i];
                double r   = length(d);
                if (r <= 0) continue;
                double phi;
                if (r < p.rHardCore)
                    phi = coreForce - p.coreStrength * (p.rHardCore - r) / p.rHardCore;
                else if (r < p.rAttract)
                    phi = 0.25 * (r - p.rEq) / (p.rAttract - p.rHardCore);
                else
                    phi = p.farAttract;
                coh += d * (phi / r);
            }
            coh *= p.cohesion / directed[i].size();
        }

        // speed control V(w)
        double wi = w[i];
        if (p.speedControl == Params::SpeedControl::Linear)
            Fwi -= p.g * (wi - 1.0);
        else {
            double q = wi * wi - 1.0;
            Fwi -= 8.0 * p.lambda * wi * q * q * q;
        }

        F[i]  = Fi + coh + h[i];
        Fw[i] = Fwi + entropic + p.cohesionSpeedGain * dot(coh, u[i]) + b[i];
    }
}

void Flock::stepOverdamped()
{
    const Params &p = mParams;
    computeForces(mF, mFw);
    const double noise = std::sqrt(2.0 * p.T * p.dt / p.eta);
    for (int i = 0; i < size(); i++) {
        Vec3 du = perp(mF[i], u[i]) * (p.dt / p.eta) + perp(gaussVec(), u[i]) * noise;
        u[i]    = normalize(u[i] + du);
    }
    stepSpeedAndPosition(mFw);
}

void Flock::stepInertial()
{
    // BAOAB splitting for the heading/spin pair: half kick, half rotation, exact
    // friction+noise (Ornstein-Uhlenbeck) on the spin, half rotation, half kick.
    const Params &p  = mParams;
    const double  hd = 0.5 * p.dt;
    auto rotateHalf = [&](int i) {
        double s = length(sigma[i]);
        if (s > 0) u[i] = normalize(rotate(u[i], sigma[i] / s, s / p.chi * hd));
    };

    computeForces(mF, mFw);
    const double c     = std::exp(-p.eta / p.chi * p.dt);
    const double noise = std::sqrt(p.chi * p.T * (1.0 - c * c));
    for (int i = 0; i < size(); i++) {
        sigma[i] += cross(u[i], mF[i]) * hd;
        rotateHalf(i);
        sigma[i] = sigma[i] * c + perp(gaussVec(), u[i]) * noise;
        rotateHalf(i);
    }

    std::vector<double> Fw = mFw;   // speed uses forces from the start of the step
    computeForces(mF, mFw);
    for (int i = 0; i < size(); i++) {
        sigma[i] += cross(u[i], mF[i]) * hd;
        sigma[i] = perp(sigma[i], u[i]);
    }
    stepSpeedAndPosition(Fw);
}

void Flock::stepSpeedAndPosition(const std::vector<double> &Fw)
{
    const Params &p     = mParams;
    const double  noise = std::sqrt(2.0 * p.T * p.dt / p.gammaS);
    for (int i = 0; i < size(); i++) {
        double wi = w[i] + Fw[i] * (p.dt / p.gammaS) + gauss() * noise;
        w[i]      = std::clamp(wi, 0.2, 3.0);
        x[i] += u[i] * (p.dt * p.v0 * w[i]);
    }
}

void Flock::step()
{
    if (mParams.turning == Params::Turning::Inertial)
        stepInertial();
    else
        stepOverdamped();
    mTime += mParams.dt;
    mSinceNeighbors += mParams.dt;
}

void Flock::advance(double seconds, double neighborInterval)
{
    mUprev = u;
    int steps = std::max(1, (int) std::lround(seconds / mParams.dt));
    for (int s = 0; s < steps; s++) {
        if (mSinceNeighbors >= neighborInterval) updateNeighbors();
        step();
    }
    double elapsed = steps * mParams.dt;
    for (int i = 0; i < size(); i++)
        mTurnRate[i] = std::acos(std::clamp(dot(mUprev[i], u[i]), -1.0, 1.0)) / elapsed;
}

} // namespace sim
