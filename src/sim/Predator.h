#pragma once
#include "Flock.h"
#include "Measure.h"
#include <random>
#include <vector>

namespace sim {

// Predator model (docs/SIMULATION_DESIGN.md section 4).
//
// The predator enters the flock's dynamics only as external fields on the birds that detect it:
// a heading field h_i pointing away from the predator and a speed field b_i (an escape speed-up).
// In the energy these are the -h.u and -b.s terms, the same role boundary birds play in the
// paper's Appendix D.2. Everything else about the response comes from the flock's own couplings.
struct PredatorParams {
    double detectRadius      = 5.0;    // only birds this close sense the predator directly (m) [ours]
    double reactionDelay     = 0.07;   // s (reported 50-90 ms; Hemelrijk 2015, Papadopoulou 2026)
    double headingField      = 1.0;    // strength, in units of J * n_c
    double speedField        = 0.3;    // strength, in units of J * n_c
    double attackSpeedFactor = 1.6;    // predator speed / v0 (reported 1.3-2x)
    double maxTurnRate       = 3.0;    // predator steering limit (rad/s)
    double startDistance     = 45.0;   // where an attack starts, from the flock center (m)
    double captureRadius     = 0.5;    // m
    double strikeDistance    = 1.5;    // after closing to this distance from its target bird, the
                                       // predator breaks away outward (an edge strike) [ours]
    double strikeTime        = 0.6;    // ...or this long after the first bird detects it [ours]
    bool   autoAttack        = false;
    double autoInterval      = 15.0;   // s between automatic attacks
    double falseAlarmRate    = 0.0;    // spontaneous single-bird startles per second (0 = off)
    double falseAlarmLength  = 0.4;    // s
};

class Predator {
  public:
    enum class State { Idle, Attacking, Passing, Retreating };

    explicit Predator(uint64_t seed = 7) : mRng(seed) {}

    PredatorParams params;

    // Starts an attack run from a random direction (above and to the side of the flock).
    void triggerAttack(const Flock &flock);
    // Moves the predator by dt and writes the fields h and b into the flock. Call once per frame
    // before flock.advance(dt).
    void update(Flock &flock, double dt);

    State state() const { return mState; }
    bool  active() const { return mState != State::Idle; }
    Vec3  position() const { return mPos; }
    Vec3  velocity() const { return mVel; }
    bool  detecting(int i) const { return i < (int) mDetectOn.size() && mDetectOn[i]; }
    bool  falseAlarm(int i) const { return i < (int) mAlarmUntil.size() && mAlarmUntil[i] > mTime; }

    // Outcome of the most recent attack (model predictions, not claims about real starlings)
    struct Report {
        int    attacks        = 0;
        double startTime      = 0;
        int    detectedEver   = 0;     // birds that sensed the predator directly
        int    captures       = 0;
        double fractionTurned = 0;     // birds that turned > 30 deg within 2 s of first detection
        double frontSpeed     = 0;     // m/s, from the TurnTracker
        double minPolarization = 1;
        int    maxFragments   = 1;
        bool   finished       = false;
    };
    const Report &report() const { return mReport; }

  private:
    void resize(int N);

    std::mt19937_64 mRng;
    State           mState = State::Idle;
    Vec3            mPos, mVel;
    int             mTarget = -1;
    double          mTime = 0, mStateTime = 0, mNextAuto = 0, mFirstDetect = -1;
    std::vector<double> mDetectSince;   // time at which bird i entered the radius (-1 = outside)
    std::vector<char>   mDetectOn, mEverDetected, mCaught;
    std::vector<double> mAlarmUntil;
    std::vector<Vec3>   mAlarmDir;
    TurnTracker         mTracker;
    Report              mReport;
};

} // namespace sim
