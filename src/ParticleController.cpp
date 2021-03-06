#include "cinder/app/App.h"
#include "cinder/Rand.h"
#include "cinder/Vector.h"
#include "cinder/gl/gl.h"
#include "ParticleController.h"

using namespace ci;
using std::list;

ParticleController::ParticleController()
{
}

void ParticleController::applyForceToParticles(float zoneRadiusSqrd)
{
    for (list<Particle>::iterator p1 = mParticles.begin(); p1 != mParticles.end(); ++p1)
    {
        list<Particle>::iterator p2 = p1;
        for (++p2; p2 != mParticles.end(); ++p2)
        {
            vec3 dir = p1->mPos - p2->mPos;
            float distSqrd = length2(dir);

            // SEPARATION
            if (distSqrd <= zoneRadiusSqrd)
            {
                float F = (zoneRadiusSqrd / distSqrd - 1.0f) * 0.01f;
                normalize(dir);
                dir *= F;
                p1->mAcc += dir;
                p2->mAcc -= dir;
            }
        }
    }
}

void ParticleController::pullToCenter(const ci::vec3 &center)
{
    for (list<Particle>::iterator p = mParticles.begin(); p != mParticles.end(); ++p)
    {
        p->pullToCenter(center);
    }
}

void ParticleController::update(bool flatten)
{
    for (list<Particle>::iterator p = mParticles.begin(); p != mParticles.end(); ++p)
    {
        p->update(flatten);
    }
}

void ParticleController::draw()
{
    gl::color(ColorA(1.0f, 1.0f, 1.0f, 1.0f));
    for (list<Particle>::iterator p = mParticles.begin(); p != mParticles.end(); ++p)
    {
        p->draw();
    }

    gl::begin(GL_LINES);
    for (list<Particle>::iterator p = mParticles.begin(); p != mParticles.end(); ++p)
    {
        p->drawTail();
    }
    gl::end();
}

void ParticleController::addParticles(int amt)
{
    for (int i = 0; i < amt; i++)
    {
        vec3 pos = Rand::randVec3() * Rand::randFloat(50.0f, 250.0f);
        vec3 vel = Rand::randVec3() * 2.0f;
        mParticles.push_back(Particle(pos, vel));
    }
}
