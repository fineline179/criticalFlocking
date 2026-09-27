#pragma once
#include "Flock.h"
#include <array>
#include <complex>
#include <deque>
#include <vector>

namespace sim {

// Snapshot statistics (docs/SIMULATION_DESIGN.md section 5)
struct Stats {
    double polarization = 0;     // |mean heading|
    Vec3   meanHeading;
    Vec3   center;
    double meanSpeed    = 0;     // m/s
    double speedSpread  = 0;     // SD of individual speeds / mean speed
    double avgNeighbors = 0;     // directed topological neighbors per bird
    double nearestDist  = 0;     // mean nearest-neighbor distance r1 (m)
    double Qint         = 0;     // similarity to neighbors, paper eq. (1), with v0 = mean speed
    double L            = 0;     // flock size: max distance between two birds (m)
    int    fragments    = 1;     // connected groups of >= kMinFragment birds
    // Connected correlation functions of heading and speed fluctuations, binned in r,
    // and their zero crossings (the correlation lengths)
    double              binWidth = 0.5;
    std::vector<double> Cdir, Csp;
    double              xiDir = 0, xiSp = 0;
};

Stats measure(const Flock &flock, bool correlations = true, double binWidth = 0.5);

// Tracks a turn spreading from an origin: when each bird's heading has deviated from the
// pre-turn mean heading by more than a threshold, its distance from the origin, and its peak
// turn rate. Gives the front speed (fit of distance vs arrival time) and attenuation.
class TurnTracker {
  public:
    // If turnDirection is nonzero, a bird counts as turned once its heading has rotated by
    // thresholdDeg *toward* that direction (robust to spontaneous wobbles); otherwise once it has
    // deviated from the pre-turn mean heading by thresholdDeg in any direction.
    void start(const Flock &flock, const Vec3 &origin, double thresholdDeg = 30.0,
               const Vec3 &turnDirection = Vec3());
    void update(const Flock &flock);
    bool active() const { return mActive; }
    void stop() { mActive = false; }

    struct Result {
        double fractionTurned = 0;
        double frontSpeed     = 0;   // m/s, slope of distance vs arrival time (all birds)
        double fitR2          = 0;
        double frontSpeedBinned = 0; // m/s, from distance-binned mean arrival times
        double attenuation    = 0;   // mean peak turn rate, farthest third / nearest third
        double elapsed        = 0;
        int    turned         = 0;
        // mean arrival time in distance bins (for diagnostics), bin width = profileBin
        double              profileBin = 2.0;
        std::vector<double> profileTime;
        std::vector<int>    profileCount;
    };
    Result result() const;

  private:
    bool                mActive = false;
    double              mT0 = 0, mCosThreshold = 0, mLastTime = 0;
    Vec3                mRefHeading, mOrigin, mTurnDir;
    double              mSinThreshold = 0;
    std::vector<double> mDist, mArrival, mPeakRate;
};

// Probes whether spontaneous heading fluctuations are overdamped or oscillatory, following the
// observable of Cavagna et al. PRL 2026: the time autocorrelation of Fourier modes of the heading
// fluctuations at wavenumbers of a few inverse meters, in the flock's center-of-mass frame.
class FluctuationProbe {
  public:
    // sampleInterval: time between calls to sample(); maxLag: longest lag analysed;
    // window: how much history to keep for averaging over time origins
    FluctuationProbe(double sampleInterval, std::vector<double> kValues = { 1.5, 3.0 },
                     double maxLag = 2.0, double window = 60.0);
    void sample(const Flock &flock);
    void clear();

    struct Result {
        std::vector<double> k;
        std::vector<double> minNormalized;   // most negative value of C(t)/C(0) per k
        std::vector<double> halfLife;        // time for C(t)/C(0) to drop below 1/2 (s)
        bool                valid = false;
    };
    Result result() const;

  private:
    double                                         mInterval;
    std::vector<double>                            mK;
    std::vector<Vec3>                              mDirs;
    size_t                                         mMaxLagSamples, mWindowSamples;
    // per k-vector, time series of the complex mode amplitude (3 components)
    std::vector<std::deque<std::array<std::complex<double>, 3>>> mSeries;
};

} // namespace sim
