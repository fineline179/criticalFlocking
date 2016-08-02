#pragma once
#include "cinder/Vector.h"
#include "cinder/Color.h"
#include <vector>

class Particle {
 public:
	Particle();
	Particle( ci::Vec3f pos, ci::Vec3f vel );
	void pullToCenter( const ci::Vec3f &center );
	void update( bool flatten );
	void limitSpeed();
    void draw(float radScale = 1.0f);
	void drawTail(float dim = 1.0f);
    inline void updatePos_spp(double x, double y, double z)
    {
        mPos.set(x, y, z); 
    };
    inline void updateVel_spp(double vx, double vy, double vz) 
    {
        mVel.set(vx, vy, vz);
    };
	
	ci::Vec3f	mPos;
	ci::Vec3f	mTailPos;
	ci::Vec3f	mVel;
	ci::Vec3f	mVelNormal;
	ci::Vec3f	mAcc;
	
	float		mDecay;
	float		mRadius;
	float		mLength;
	float		mMaxSpeed, mMaxSpeedSqrd;
	float		mMinSpeed, mMinSpeedSqrd;

	bool		mIsDead;
};