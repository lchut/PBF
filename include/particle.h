#ifndef PBF_PARTICLE_H_
#define PBF_PARTICLE_H_

#include "config.h"
#include <vector>
#include <cmath>
#include "AABB.h"

struct Particle {
    Particle() : velocity(0.0f) {}
    Particle(const glm::vec3& v) : velocity(v) {}
    Particle(const glm::vec3& pos, const glm::vec3& v) : position(pos), velocity(v) {}
    glm::vec3 position;
    glm::vec3 velocity;
};

class ParticleManager {
public:
    ParticleManager() = default;
    ParticleManager(uint32_t nParticels) {
        mParticles = new Particle[nParticels];
    }
    Particle* getParticlesPtr(){ return mParticles; }

    Particle*  mParticles;

    float radius;
};

#endif