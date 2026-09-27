// Headless driver for the simulation core, for calibration and validation runs.
//
//   flocksim_cli <stats|turn|speedpush|attack|probe> [--turning overdamped|inertial|nonlinear]
//                [--speed stiff|critical|marginal] [--N 300] [--seed 1] [--warmup 20]
//                [--measure 20] [--set name=value ...]
//
//   stats: time-averaged snapshot statistics (docs/SIMULATION_DESIGN.md section 5)
//   turn:  a few birds at the flock's edge are pushed into a sustained 90-degree turn; reports
//          how the turn front spreads
//   probe: whether spontaneous heading fluctuations are overdamped or oscillatory
#include "sim/Flock.h"
#include "sim/Measure.h"
#include "sim/Predator.h"
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <string>

using namespace sim;

namespace {

struct Options {
    std::string   command  = "stats";
    TurningPreset turning  = TurningPreset::Inertial;
    SpeedPreset   speed    = SpeedPreset::Stiff;
    int           N        = 300;
    uint64_t      seed     = 1;
    double        warmup   = 20;
    double        duration = 20;
    // turn test
    int    initiators  = 3;
    double fieldFactor = 2.0;   // initiator field strength in units of J * n_c
    double turnTime    = 3.0;
    double pulse       = 0.0;   // if > 0, the initiator field is only on for this long
    std::map<std::string, double> overrides;
};

bool setParam(Params &p, const std::string &name, double v)
{
    std::map<std::string, double *> fields = {
        { "v0", &p.v0 },       { "chi", &p.chi },        { "eta", &p.eta },
        { "J", &p.J },         { "J4", &p.J4 },          { "T", &p.T },
        { "g", &p.g },         { "lambda", &p.lambda },  { "gammaS", &p.gammaS },
        { "speedCoupling", &p.speedCoupling },           { "balanceAngleDeg", &p.balanceAngleDeg },
        { "rHardCore", &p.rHardCore }, { "rEq", &p.rEq }, { "rAttract", &p.rAttract },
        { "cohesion", &p.cohesion },   { "dt", &p.dt },
        { "farAttract", &p.farAttract }, { "coreStrength", &p.coreStrength },
        { "cohesionSpeedGain", &p.cohesionSpeedGain },
    };
    if (name == "nc") { p.nc = (int) v; return true; }
    if (name == "nearest") {
        p.neighborRule = v != 0 ? Params::NeighborRule::Nearest : Params::NeighborRule::Balanced;
        return true;
    }
    auto it = fields.find(name);
    if (it == fields.end()) return false;
    *it->second = v;
    return true;
}

void usage()
{
    std::fprintf(stderr, "usage: flocksim_cli <stats|turn|speedpush|attack|probe> [--turning overdamped|inertial|nonlinear]"
                         " [--speed stiff|critical|marginal] [--N n] [--seed s] [--warmup sec]"
                         " [--measure sec] [--initiators n] [--field f] [--turn-time sec] [--pulse sec]"
                         " [--set name=value ...]\n");
    std::exit(1);
}

Options parse(int argc, char **argv)
{
    Options o;
    if (argc < 2) usage();
    o.command = argv[1];
    for (int a = 2; a < argc; a++) {
        std::string k = argv[a];
        auto next = [&]() -> std::string { if (a + 1 >= argc) usage(); return argv[++a]; };
        if (k == "--turning") {
            std::string v = next();
            o.turning = v == "overdamped" ? TurningPreset::Overdamped
                      : v == "nonlinear"  ? TurningPreset::NonlinearInertial
                                          : TurningPreset::Inertial;
        }
        else if (k == "--speed") {
            std::string v = next();
            o.speed = v == "critical" ? SpeedPreset::NearCritical
                    : v == "marginal" ? SpeedPreset::Marginal
                                      : SpeedPreset::Stiff;
        }
        else if (k == "--N") o.N = std::atoi(next().c_str());
        else if (k == "--seed") o.seed = std::strtoull(next().c_str(), nullptr, 10);
        else if (k == "--warmup") o.warmup = std::atof(next().c_str());
        else if (k == "--measure") o.duration = std::atof(next().c_str());
        else if (k == "--initiators") o.initiators = std::atoi(next().c_str());
        else if (k == "--field") o.fieldFactor = std::atof(next().c_str());
        else if (k == "--turn-time") o.turnTime = std::atof(next().c_str());
        else if (k == "--pulse") o.pulse = std::atof(next().c_str());
        else if (k == "--set") {
            std::string kv = next();
            auto        eq = kv.find('=');
            if (eq == std::string::npos) usage();
            o.overrides[kv.substr(0, eq)] = std::atof(kv.c_str() + eq + 1);
        }
        else usage();
    }
    return o;
}

void printParams(const Params &p)
{
    std::printf("params: turning=%s J=%g J4=%g eta=%g chi=%g T=%g speed=%s g=%g lambda=%g "
                "gammaS=%g cohesion=%g N=%d dt=%g\n",
                p.turning == Params::Turning::Inertial ? "inertial" : "overdamped", p.J, p.J4,
                p.eta, p.chi, p.T,
                p.speedControl == Params::SpeedControl::Linear ? "linear" : "marginal", p.g,
                p.lambda, p.gammaS, p.cohesion, p.N, p.dt);
}

} // namespace

int main(int argc, char **argv)
{
    Options o = parse(argc, argv);
    Params  p = makePreset(o.turning, o.speed, o.N);
    for (auto &[k, v] : o.overrides)
        if (!setParam(p, k, v)) { std::fprintf(stderr, "unknown parameter %s\n", k.c_str()); return 1; }
    printParams(p);

    Flock flock(p, o.seed);
    const double frame = 1.0 / 60.0;
    auto         wall0 = std::chrono::steady_clock::now();
    for (double t = 0; t < o.warmup; t += frame) flock.advance(frame);

    if (o.command == "stats") {
        double n = 0, phi = 0, spd = 0, spread = 0, nb = 0, L = 0, xiD = 0, xiS = 0, frag = 0, r1 = 0,
               qint = 0;
        double phiMin = 1;
        for (double t = 0; t < o.duration; t += frame) {
            flock.advance(frame);
            if (std::fmod(t, 0.5) < frame) {
                Stats s = measure(flock);
                phi += s.polarization; spd += s.meanSpeed; spread += s.speedSpread;
                nb += s.avgNeighbors; L += s.L; xiD += s.xiDir; xiS += s.xiSp;
                frag += s.fragments; phiMin = std::min(phiMin, s.polarization);
                r1 += s.nearestDist;
                qint += s.Qint;
                n++;
            }
        }
        std::printf("polarization=%.4f polarization_min=%.4f mean_speed=%.2f speed_spread=%.3f "
                    "neighbors=%.2f L=%.2f xi_dir=%.2f xi_sp=%.2f xi_dir_over_L=%.3f "
                    "xi_sp_over_L=%.3f fragments=%.2f r1=%.2f Qint=%.4f\n",
                    phi / n, phiMin, spd / n, spread / n, nb / n, L / n, xiD / n, xiS / n,
                    xiD / L, xiS / L, frag / n, r1 / n, qint / n);
    }
    else if (o.command == "turn") {
        Stats s0 = measure(flock, false);
        // lateral edge: the birds farthest along an axis perpendicular to the mean heading
        Vec3 side = normalize(perp(Vec3(0, 1, 0), s0.meanHeading));
        std::vector<int> order(flock.size());
        for (int i = 0; i < flock.size(); i++) order[i] = i;
        std::sort(order.begin(), order.end(), [&](int a, int b) {
            return dot(flock.x[a] - s0.center, side) > dot(flock.x[b] - s0.center, side);
        });
        Vec3 origin;
        for (int k = 0; k < o.initiators; k++) origin += flock.x[order[k]];
        origin = origin / o.initiators;

        // sustained field turning the initiators toward the outward side (a 90-degree turn)
        double      hmag = o.fieldFactor * p.J * p.nc;
        TurnTracker tracker;
        tracker.start(flock, origin, 30.0, side);
        double fragMax = 1;
        for (double t = 0; t < o.turnTime; t += frame) {
            bool on = o.pulse <= 0 || t < o.pulse;
            for (int k = 0; k < o.initiators; k++) flock.h[order[k]] = on ? side * hmag : Vec3();
            flock.advance(frame);
            tracker.update(flock);
            fragMax = std::max(fragMax, (double) measure(flock, false).fragments);
        }
        auto r = tracker.result();
        Stats s1 = measure(flock, false);
        std::printf("fraction_turned=%.3f front_speed=%.2f front_speed_binned=%.2f "
                    "binned_over_v0=%.2f fit_r2=%.3f attenuation=%.3f fragments_max=%.0f "
                    "polarization_after=%.3f L=%.2f\n",
                    r.fractionTurned, r.frontSpeed, r.frontSpeedBinned, r.frontSpeedBinned / p.v0,
                    r.fitR2, r.attenuation, fragMax, s1.polarization, s0.L);
        std::printf("arrival profile (distance_bin_start: mean_arrival_s n):");
        for (size_t b = 0; b < r.profileTime.size(); b++)
            if (r.profileCount[b])
                std::printf(" %.0fm:%.3f(%d)", b * r.profileBin, r.profileTime[b], r.profileCount[b]);
        std::printf("\n");
    }
    else if (o.command == "speedpush") {
        // Sustained speed-up field on a few edge birds; how far does the speed change spread?
        Stats s0   = measure(flock, false);
        Vec3  side = normalize(perp(Vec3(0, 1, 0), s0.meanHeading));
        std::vector<int> order(flock.size());
        for (int i = 0; i < flock.size(); i++) order[i] = i;
        std::sort(order.begin(), order.end(), [&](int a, int b) {
            return dot(flock.x[a] - s0.center, side) > dot(flock.x[b] - s0.center, side);
        });
        // baseline: each bird's mean speed over 2 s; response: over the last 2 s of a 4 s push
        auto meanSpeeds = [&](double seconds, bool push) {
            std::vector<double> m(flock.size(), 0.0);
            int n = 0;
            for (double t = 0; t < seconds; t += frame, n++) {
                for (int k = 0; k < o.initiators; k++)
                    flock.b[order[k]] = push ? o.fieldFactor * p.J * p.nc : 0.0;
                flock.advance(frame);
                for (int i = 0; i < flock.size(); i++) m[i] += flock.w[i];
            }
            for (double &v : m) v /= n;
            return m;
        };
        std::vector<double> base = meanSpeeds(2.0, false);
        Vec3 origin;
        for (int k = 0; k < o.initiators; k++) origin += flock.x[order[k]];
        origin = origin / o.initiators;
        std::vector<double> d(flock.size());
        for (int i = 0; i < flock.size(); i++) d[i] = length(flock.x[i] - origin);
        meanSpeeds(2.0, true);
        std::vector<double> resp = meanSpeeds(2.0, true);
        const double bin = 3.0;
        std::vector<double> dv(40, 0.0);
        std::vector<int> cnt(40, 0);
        double flockMean = 0;
        for (int i = 0; i < flock.size(); i++) {
            double dw = (resp[i] - base[i]) * p.v0;
            flockMean += dw / flock.size();
            int b = std::min(39, (int) (d[i] / bin));
            dv[b] += dw;
            cnt[b]++;
        }
        std::printf("speed change (m/s) by distance from pushed birds:");
        for (int b = 0; b < 40; b++)
            if (cnt[b] >= 3) std::printf(" %.0fm:%+.2f(%d)", b * bin, dv[b] / cnt[b], cnt[b]);
        std::printf("\nflock_mean_speed_change=%.2f L=%.2f\n", flockMean, s0.L);
    }
    else if (o.command == "attack") {
        // One predator attack run; reports the outcome (model predictions)
        Predator pred(o.seed + 100);
        pred.params.headingField = o.fieldFactor;
        pred.triggerAttack(flock);
        for (double t = 0; t < 12.0 && !pred.report().finished; t += frame) {
            pred.update(flock, frame);
            flock.advance(frame);
        }
        const auto &r = pred.report();
        std::printf("detected_directly=%d detected_fraction=%.3f fraction_turned_2s=%.3f front_speed=%.2f captures=%d "
                    "min_polarization=%.3f max_fragments=%d finished=%d\n",
                    r.detectedEver, double(r.detectedEver) / flock.size(), r.fractionTurned,
                    r.frontSpeed, r.captures, r.minPolarization,
                    r.maxFragments, (int) r.finished);
    }
    else if (o.command == "probe") {
        const double     interval = 0.02;
        FluctuationProbe probe(interval, { 0.75, 1.5, 3.0 }, 3.0, o.duration);
        double           since = 0;
        for (double t = 0; t < o.duration; t += frame) {
            flock.advance(frame);
            since += frame;
            if (since >= interval) { probe.sample(flock); since = 0; }
        }
        auto r = probe.result();
        for (size_t a = 0; a < r.k.size(); a++)
            std::printf("k=%.2f min_normalized_autocorr=%.3f half_life=%.3f %s\n", r.k[a],
                        r.minNormalized[a], r.halfLife[a],
                        r.minNormalized[a] < -0.1 ? "OSCILLATORY" : "overdamped");
    }
    else usage();

    double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - wall0).count();
    std::printf("sim_time=%.1f wall_time=%.2f realtime_factor=%.1f\n", flock.time(), wall,
                flock.time() / wall);
    return 0;
}
