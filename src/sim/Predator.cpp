#include "Predator.h"
#include <algorithm>
#include <cmath>

namespace sim {

namespace {
const Vec3 kUp(0, 1, 0);   // world up (matches the app's camera)

// turn v toward desired direction by at most maxAngle, keeping |v|
Vec3 steer(const Vec3 &v, const Vec3 &desired, double maxAngle, double speed)
{
    Vec3   a     = normalize(v), d = normalize(desired);
    double angle = std::acos(std::clamp(dot(a, d), -1.0, 1.0));
    if (angle <= maxAngle) return d * speed;
    Vec3 axis = cross(a, d);
    if (length2(axis) < 1e-12) axis = normalize(perp(kUp, a));
    return normalize(rotate(a, normalize(axis), maxAngle)) * speed;
}
} // namespace

void Predator::resize(int N)
{
    if ((int) mDetectSince.size() == N) return;
    mDetectSince.assign(N, -1.0);
    mDetectOn.assign(N, 0);
    mEverDetected.assign(N, 0);
    mCaught.assign(N, 0);
    mAlarmUntil.assign(N, -1.0);
    mAlarmDir.assign(N, Vec3());
}

void Predator::triggerAttack(const Flock &flock)
{
    resize(flock.size());
    Stats s = measure(flock, false);

    // approach from a random direction around the flock's heading, biased from above
    std::uniform_real_distribution<double> ang(0.0, 2.0 * 3.14159265358979);
    Vec3   side1 = normalize(perp(kUp, s.meanHeading));
    Vec3   side2 = normalize(cross(s.meanHeading, side1));
    double a     = ang(mRng);
    Vec3   dir   = normalize(side1 * std::cos(a) + side2 * std::sin(a) + kUp * 0.6);
    mPos         = s.center + dir * params.startDistance + s.meanHeading * 10.0;

    // target the bird closest to the predator (an edge bird on the approach side)
    double best = 1e300;
    for (int i = 0; i < flock.size(); i++) {
        double d2 = length2(flock.x[i] - mPos);
        if (d2 < best) { best = d2; mTarget = i; }
    }
    mVel        = normalize(flock.x[mTarget] - mPos) * (params.attackSpeedFactor * flock.params().v0);
    mState      = State::Attacking;
    mStateTime  = 0;
    mFirstDetect = -1;
    std::fill(mEverDetected.begin(), mEverDetected.end(), 0);
    std::fill(mCaught.begin(), mCaught.end(), 0);
    int attacks        = mReport.attacks + 1;
    mReport            = Report();
    mReport.attacks    = attacks;
    mReport.startTime  = mTime;
    mReport.minPolarization = s.polarization;
    mTracker.stop();
}

void Predator::update(Flock &flock, double dt)
{
    const int     N = flock.size();
    const Params &fp = flock.params();
    resize(N);
    mTime += dt;
    mStateTime += dt;

    if (params.autoAttack && mState == State::Idle && mTime >= mNextAuto) {
        triggerAttack(flock);
        mNextAuto = mTime + params.autoInterval;
    }

    // predator motion
    const double speed = params.attackSpeedFactor * fp.v0;
    if (mState != State::Idle) {
        Stats s = measure(flock, false);
        switch (mState) {
            case State::Attacking: {
                Vec3 toTarget = flock.x[mTarget] - mPos;
                mVel          = steer(mVel, toTarget, params.maxTurnRate * dt, speed);
                // strike at the edge bird, then break away (don't plow through the flock)
                bool strikeOver = mFirstDetect >= 0 && mTime - mFirstDetect > params.strikeTime;
                if (length(toTarget) < params.strikeDistance || strikeOver || mStateTime > 10.0) {
                    mState     = State::Passing;
                    mStateTime = 0;
                }
                break;
            }
            case State::Passing: {
                Vec3 outward = normalize(mPos - s.center);
                mVel         = steer(mVel, outward + mVel / speed, params.maxTurnRate * dt, speed);
                if (mStateTime > 1.5) { mState = State::Retreating; mStateTime = 0; }
                break;
            }
            case State::Retreating: {
                Vec3 away = normalize(mPos - s.center) + kUp;
                mVel      = steer(mVel, away, params.maxTurnRate * dt, speed);
                if (mStateTime > 4.0 || length(mPos - s.center) > 80.0) {
                    mState           = State::Idle;
                    mReport.finished = true;
                }
                break;
            }
            case State::Idle: break;
        }
        mPos += mVel * dt;

        if (mState != State::Idle) {
            mReport.minPolarization = std::min(mReport.minPolarization, s.polarization);
            mReport.maxFragments    = std::max(mReport.maxFragments, s.fragments);
        }
    }

    // detection (with reaction delay) and the resulting fields
    const double hMag = params.headingField * fp.J * fp.nc;
    const double bMag = params.speedField * fp.J * fp.nc;
    for (int i = 0; i < N; i++) {
        flock.h[i] = Vec3();
        flock.b[i] = 0.0;
        bool inRange = mState != State::Idle && length(flock.x[i] - mPos) < params.detectRadius;
        if (!inRange)
            mDetectSince[i] = -1;
        else if (mDetectSince[i] < 0)
            mDetectSince[i] = mTime;
        mDetectOn[i] = inRange && mTime - mDetectSince[i] >= params.reactionDelay;

        if (mDetectOn[i]) {
            flock.h[i] = normalize(flock.x[i] - mPos) * hMag;
            flock.b[i] = bMag;
            if (!mEverDetected[i]) {
                mEverDetected[i] = 1;
                mReport.detectedEver++;
            }
            if (mFirstDetect < 0) {
                mFirstDetect = mTime;
                mTracker.start(flock, flock.x[i]);
            }
        }
        if ((mState == State::Attacking || mState == State::Passing) && !mCaught[i] &&
            length(flock.x[i] - mPos) < params.captureRadius) {
            mCaught[i] = 1;
            mReport.captures++;
        }
    }

    // spontaneous false alarms: a random bird startles as if it had seen something
    if (params.falseAlarmRate > 0) {
        std::uniform_real_distribution<double> uni(0.0, 1.0);
        if (uni(mRng) < params.falseAlarmRate * dt) {
            int i = std::min(N - 1, (int) (uni(mRng) * N));
            std::normal_distribution<double> g;
            mAlarmUntil[i] = mTime + params.falseAlarmLength;
            mAlarmDir[i]   = normalize(perp(Vec3(g(mRng), g(mRng), g(mRng)), flock.u[i]));
        }
        for (int i = 0; i < N; i++)
            if (mAlarmUntil[i] > mTime) {
                flock.h[i] += mAlarmDir[i] * hMag;
                flock.b[i] += bMag;
            }
    }

    // turn-front measurement, from the first detection
    if (mTracker.active()) {
        mTracker.update(flock);
        auto r = mTracker.result();
        if (r.elapsed <= 2.0) mReport.fractionTurned = r.fractionTurned;
        mReport.frontSpeed = r.frontSpeedBinned;
        if (r.elapsed > 4.0) mTracker.stop();
    }
}

} // namespace sim
