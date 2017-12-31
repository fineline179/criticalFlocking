#pragma once
#include "cinder/Vector.h"
#include "cinder/Color.h"
#include <vector>

class Particle {
public:
    Particle();
    Particle(ci::vec3 pos, ci::vec3 vel);
    void pullToCenter(const ci::vec3 &center);
    void update(bool flatten);
    void limitSpeed();
    void draw(float radScale = 1.0f);
    void drawTail(float dim = 1.0f);
    inline void updatePos_spp(double x, double y, double z)
    {
        mPos = ci::vec3(x, y, z);
    };
    inline void updateVel_spp(double vx, double vy, double vz)
    {
        mVel = ci::vec3(vx, vy, vz);
    };

    ci::vec3	mPos;
    ci::vec3	mTailPos;
    ci::vec3	mVel;
    ci::vec3	mVelNormal;
    ci::vec3	mAcc;

    float		mDecay;
    float		mRadius;
    float		mLength;
    ci::vec3    mCubeSize;
    float		mMaxSpeed, mMaxSpeedSqrd;
    float		mMinSpeed, mMinSpeedSqrd;

    bool		mIsDead;
};