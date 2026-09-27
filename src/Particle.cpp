#include "Particle.h"
#include "cinder/Rand.h"
#include "cinder/gl/gl.h"
#include "cinder/app/App.h"

using namespace ci;

Particle::Particle()
{
}

Particle::Particle(vec3 pos, vec3 vel)
{
	mPos			= pos;
	mTailPos		= pos;
	mVel			= vel;
	mVelNormal		= vec3(0.0, 1.0, 0.0);
	mAcc			= vec3(0.0, 0.0, 0.0);
	
	mMaxSpeed		= Rand::randFloat( 2.0f, 3.0f );
	mMaxSpeedSqrd	= mMaxSpeed * mMaxSpeed;
	mMinSpeed		= Rand::randFloat( 1.0f, 1.5f );
	mMinSpeedSqrd	= mMinSpeed * mMinSpeed;
	
	mDecay			= 0.99f;
    mRadius         = 0.05f;/*2.0f;*/
	mLength			= 0.25f /*10.0f*/;

}

void Particle::pullToCenter(const vec3 &center)
{
	vec3 dirToCenter = mPos - center;
    float distToCenter = dirToCenter.length();
	float maxDistance = 300.0f;
	
    if (distToCenter > maxDistance)
    {
        dirToCenter = normalize(dirToCenter);
		float pullStrength = 0.0001f;
        mVel -= dirToCenter * ((distToCenter - maxDistance) * pullStrength);
	}
}	

void Particle::update(bool flatten)
{	
    if (flatten) mAcc.z = 0.0f;
    mVel += mAcc;
    mVelNormal = normalize(mVel);
    limitSpeed();
	
    mPos += mVel;
    mTailPos = mPos - mVelNormal * mLength;
	
    if (flatten) mPos.z = 0.0f;
		
    mVel *= mDecay;
    mAcc = vec3(0.0, 0.0, 0.0);
}

void Particle::limitSpeed()
{
    float vLengthSqrd = length2(mVel);

    if (vLengthSqrd > mMaxSpeedSqrd)
    {
        mVel = mVelNormal * mMaxSpeed;
    }
    else if (vLengthSqrd < mMinSpeedSqrd)
    {
        mVel = mVelNormal * mMinSpeed;
	}
}

void Particle::draw(const gl::BatchRef &sphereBatch, float radScale)
{
    // sphereBatch is a unit sphere. Drawing it with a model matrix is equivalent to
    // gl::drawSphere(mPos, radScale * mRadius, 8), but avoids drawSphere's per-call vertex
    // upload, which stalls the GPU on macOS and dropped the frame rate to ~6 fps with a few
    // hundred birds.
    gl::ScopedModelMatrix scpModel;
    gl::translate(mPos);
    gl::scale(vec3(radScale * mRadius));
    sphereBatch->draw();
}

void Particle::drawTail(float dim)
{
    gl::color(ColorA(dim, dim, dim, 1.0f));
    gl::vertex(mPos);
    gl::color(ColorA(dim, 0.0f, 0.0f, 1.0f));
    gl::vertex(mTailPos);
}
