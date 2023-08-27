#ifndef PBF_SOLVER_H_
#define PBF_SOLVER_H_

#include <omp.h>
#include <cmath>
#include <functional>
#include <mutex>
#include <atomic>
#include <memory>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <chrono>
#include <thread>

#include "particle.h"
#include "AABB.h"
#include "solverBase.h"

const float invPi = 0.31830988618;
const float delta = 0.02;
const uint32_t maxNeighborCnt = 63;

inline float equalZero(const float v) {
    return (v > -1e-3 && v < 1e-3);
}
inline float distance(const glm::vec3& v) {
    return glm::length(v);
}
inline float distance(const glm::vec3& p1, const glm::vec3& p2) {
    return glm::length(p1 - p2);
}

inline float distSquare(const glm::vec3& v) {
    return v.x * v.x + v.y * v.y + v.z * v.z;
}
inline float distSquare(const glm::vec3& p1, const glm::vec3& p2) {
    glm::vec3 v = p1 - p2;
    return v.x * v.x + v.y * v.y + v.z * v.z;
}

inline float clamp(float L, float R, float val) {
    if (val < L) { return L; }
    else if (val > R) { return R;}
    return val;
}

inline float Wpoly6(const glm::vec3& r, const float h) {
    if (glm::length(r) > h) { return 0; }
    float c = 315.0f * invPi / (64.0f * std::pow(h, 9));
    return c * std::pow(h * h - glm::dot(r, r), 3);
}

inline glm::vec3 delWspiky(const glm::vec3& r, const float h) {
    float rLen = glm::length(r);
    if (rLen > h || rLen < 0.01f) { return glm::vec3(0.0); }
    float c = -45.0f * invPi / std::pow(h, 6);
    return c * float(std::pow((h - glm::length(r)), 2)) * glm::normalize(r);
}

inline float frac(float v) {
    return v - (int)v;
}

inline float randomOffset(const glm::vec2& uv) {
    float noise = frac(std::sin(glm::dot(uv, glm::vec2(12.9898, 78.233) * 2.0f)) * 43758.5453);
    return std::abs(noise);
}

inline void print(const glm::vec3& v) {
    printf("%f %f %f\n", v.x, v.y, v.z);
}

struct Grid {
    Grid() : startIdx(0), endIdx(0) {}
    uint32_t startIdx;
    uint32_t endIdx;
};



class PBFSolver {
public:
    PBFSolver() = default;
    PBFSolver(const SolverParameter& _para, const AABB& box,
        const std::function<glm::vec3(const glm::vec3&)>& _fext =
        [](const glm::vec3& x) { return glm::vec3(0.0, 0.0, 0.0); }) :
        para(_para), fext(_fext) {
            lambda.resize(para.nParticles, 0.0f);
            particlesRho.resize(para.nParticles, 0.0f);
            particlesOffset.resize(para.nParticles, glm::vec3(0.0));
            particlesCurl.resize(para.nParticles, glm::vec3(0.0));
            particlesViscosity.resize(para.nParticles, glm::vec3(0.0));
            particlesVorticity.resize(para.nParticles, glm::vec3(0.0));
            oldPosition.resize(para.nParticles, glm::vec3(0.0));
            boundaryBox = AABB(box.pMin + glm::vec3(para.radius), box.pMax - glm::vec3(para.radius));

            gridsCnt = gridDimSize[0] * gridDimSize[1] * gridDimSize[2];
            grids = new Grid[gridsCnt];
            gridParticlesCnt = new int[gridsCnt];
            memset(gridParticlesCnt, 0, sizeof(gridParticlesCnt));
            particlesGid = new int[para.nParticles];
            gridParticlesID = new int[para.nParticles];
            neighbors = new int*[para.nParticles];
            for (int i = 0; i < para.nParticles; ++i) {
                neighbors[i] = new int[maxNeighborCnt + 1];
            }
            pMin = box.pMin;
        }
    void solve(Particle* particles);
private:
    inline void BoundParticle(Particle& particle) {
        //printf("box:");
        //print(boundaryBox.pMin);
        //print(boundaryBox.pMax);
        //printf("pos:");
        //print(particle.position);
        if (particle.position.x < boundaryBox.pMin.x) {
            particle.position.x = boundaryBox.pMin.x + delta * randomOffset(glm::vec2(particle.position.y, particle.position.z));
        }
        if (particle.position.x > boundaryBox.pMax.x) {
            particle.position.x = boundaryBox.pMax.x - delta * randomOffset(glm::vec2(particle.position.y, particle.position.z));
        }
        if (particle.position.y < boundaryBox.pMin.y) {
            particle.position.y = boundaryBox.pMin.y + delta * randomOffset(glm::vec2(particle.position.x, particle.position.z));
        }
        if (particle.position.y > boundaryBox.pMax.y) {
            particle.position.y = boundaryBox.pMax.y - delta * randomOffset(glm::vec2(particle.position.x, particle.position.z));
        }
        if (particle.position.z < boundaryBox.pMin.z) {
            particle.position.z = boundaryBox.pMin.z + delta * randomOffset(glm::vec2(particle.position.x, particle.position.y));
        }
        if (particle.position.z > boundaryBox.pMax.z) {
            particle.position.z = boundaryBox.pMax.z - delta * randomOffset(glm::vec2(particle.position.x, particle.position.y));
        }
    }
    inline float Scorr(float wval) const {
        return -para.k * std::pow(wval / Wpoly6(glm::vec3(para.deltaq, 0.0, 0.0), para.h), para.n);
    }
    inline glm::vec3 fVorticity(const glm::vec3& N, const glm::vec3& w) const {
        return 0.5f * glm::cross(N, w);
    }
    inline void collsionProcess(Particle& particle) {
        if (particle.position.x < boundaryBox.pMin.x) {
            particle.position.x = boundaryBox.pMin.x + delta * randomOffset(glm::vec2(particle.position.y, particle.position.z));
        }
        if (particle.position.x > boundaryBox.pMax.x) {
            particle.position.x = boundaryBox.pMax.x - delta * randomOffset(glm::vec2(particle.position.y, particle.position.z));
        }
        if (particle.position.y < boundaryBox.pMin.y) {
            particle.position.y = boundaryBox.pMin.y + delta * randomOffset(glm::vec2(particle.position.x, particle.position.z));
        }
        if (particle.position.y > boundaryBox.pMax.y) {
            particle.position.y = boundaryBox.pMax.y - delta * randomOffset(glm::vec2(particle.position.x, particle.position.z));
        }
        if (particle.position.z < boundaryBox.pMin.z) {
            particle.position.z = boundaryBox.pMin.z + delta * randomOffset(glm::vec2(particle.position.x, particle.position.y));
        }
        if (particle.position.z > boundaryBox.pMax.z) {
            particle.position.z = boundaryBox.pMax.z - delta * randomOffset(glm::vec2(particle.position.x, particle.position.y));
        }
    }

    inline void init() {
        for (int i = 0; i < para.nParticles; ++i) { particlesCurl[i] = glm::vec3(0.0); }
        memset(gridParticlesCnt, 0, sizeof(gridParticlesCnt));
    }

    inline bool OutOfRange(const int dim[3]) const {
        return dim[0] < 0 || dim[0] >= gridDimSize[0] ||
            dim[1] < 0 || dim[1] >= gridDimSize[1] ||
            dim[2] < 0 || dim[2] >= gridDimSize[2];
    }

    inline int computeGid(const glm::vec3& position, float invH) const {
        int dim[3];
        computeDimID(position, invH, &dim[0], &dim[1], &dim[2]);
        if (OutOfRange(dim)) {
            return -1;
        }
        else {
            return dim[0] + dim[1] * gridDimSize[0] + dim[2] * gridDimSize[0] * gridDimSize[1];
        }
    }

    inline int computeGid(const int dim[3]) const {
        if (OutOfRange(dim)) {
            return -1;
        }
        else {
            return dim[0] + dim[1] * gridDimSize[0] + dim[2] * gridDimSize[0] * gridDimSize[1];
        }
    }
    inline void computeDimID(const glm::vec3& position, float invH,  int* dimX, int* dimY, int* dimZ) const {
        *dimX = int((position.x - pMin.x) * invH);
        *dimY = int((position.y - pMin.y) * invH);
        *dimZ = int((position.z - pMin.z) * invH);
    }
    void buildNeiborhood(const Particle* particles, int N, float h);
    std::function<glm::vec3(const glm::vec3&)> fext;

    SolverParameter para;
    AABB boundaryBox;

    std::vector<glm::vec3> oldPosition;

    std::vector<float> lambda;
    std::vector<float> particlesRho;
    std::vector<glm::vec3> particlesOffset;
    std::vector<glm::vec3> particlesCurl;
    std::vector<glm::vec3> particlesVorticity;
    std::vector<glm::vec3> particlesViscosity;
    int gridDimSize[3] = {UNIFORM_GRID_DIM_SIZE_X, UNIFORM_GRID_DIM_SIZE_Y, UNIFORM_GRID_DIM_SIZE_Z};
    int gridsCnt;
    glm::vec3 pMin;
    Grid* grids;
    int* gridParticlesCnt;
    int* particlesGid;
    int* gridParticlesID;
    int** neighbors;
};


#endif