#ifndef PBF_SOLVERBASE_H_
#define PBF_SOLVERBASE_H_

#define INVPI 0.31830988618
#define RANDOM_OFFSET_COEFF 0.01
#define MAX_NEIGHBOR 63
#define BLOCK_DIM_X 32
#define BLOCK_DIM_Y 32
#define PREFIX_BLOCK_DIM_X 16
#define PREFIX_BLOCK_DIM_Y 16
#define UNIFORM_GRID_DIM_SIZE_X 64
#define UNIFORM_GRID_DIM_SIZE_Y 64
#define UNIFORM_GRID_DIM_SIZE_Z 64
#define SHARED_MEMORY_SIZE 1024

struct SolverParameter {
    // Kernel Function parameter
    float h;
    // number of particles
    int nParticles;
    // rest desnity rho0
    float restRho;
    // relaxation parameter
    float epsilon;
    // pressure coeff
    float k;
    float n;
    float deltaq;
    // Viscosity parameter
    float c;
    // delta time
    float dt;
    float invDt;
    // solver iterations
    int iterations;
    // vorticity confinement epsilon
    float vorticityConfinementEpilon;
    // particle radius
    float radius;
};

struct AABB_;
struct Vec3;
struct Grid_;
struct SolverData {
    Vec3* position;
    Vec3* velocity;
    Vec3* oldPosition;
    float* lambda;
    Vec3* particlesOffset;
    Vec3* particlesCurl;
    Vec3* particlesVorticity;
    Vec3* particlesViscosity;
    int* neighbors;
    AABB_* gridSpace;
    AABB_* boundary;

    Grid_* grid;
    int* prefixSum;
    int* blockSum;
    int* gridParticlesCnt;
    int* gridParticlesID;

};

#endif