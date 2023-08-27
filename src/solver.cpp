 #include "solver.h"

void PBFSolver::buildNeiborhood(const Particle* particles, int N, float h) {
    // calculate particle cell id
    float invH = 1.0 / h;
    #pragma parallel for 
    for (int i = 0; i < N; ++i) {
        particlesGid[i]  = computeGid(particles[i].position, invH);
    }
    
    for (int i = 0; i < N; ++i) {
        if (particlesGid[i] >= 0) {
            ++gridParticlesCnt[particlesGid[i]];
        }
    }
    int preSum = 0, tmp;
    for (int i = 0; i < gridsCnt; ++i) {
        grids[i].startIdx = preSum;
        preSum = preSum + gridParticlesCnt[i];
        grids[i].endIdx = preSum;
    }
    
    for (int i = 0; i < N; ++i) {
        int gid = particlesGid[i];
        if (gid >= 0) {
            int offset = --gridParticlesCnt[gid];
            gridParticlesID[grids[gid].startIdx + offset] = i;
        }
    }
    float h2 = h * h;
    #pragma parallel for 
    for (int i = 0; i < N; ++i) {
        int neighborCnt = 0;
        int dim[3];
        computeDimID(particles[i].position, invH, &dim[0], &dim[1], &dim[2]);
        int neighborDim[3];
        for (int xOffset = -1; xOffset <= 1 && neighborCnt < maxNeighborCnt; ++xOffset) {
            for (int yOffset = -1; yOffset <= 1 && neighborCnt < maxNeighborCnt; ++yOffset) {
                for (int zOffset = -1; zOffset <= 1 && neighborCnt < maxNeighborCnt; ++zOffset) {
                    neighborDim[0] = dim[0] + xOffset;
                    neighborDim[1] = dim[1] + yOffset;
                    neighborDim[2] = dim[2] + zOffset;
                    int gid = computeGid(neighborDim);
                    if (gid < 0) { continue; }
                    for (int j = grids[gid].startIdx; j < grids[gid].endIdx && neighborCnt < maxNeighborCnt; ++j) {
                        int nid = gridParticlesID[j];
                        if (nid != i) {
                            if (distSquare(particles[i].position - particles[nid].position) <= h2) {
                                neighbors[i][neighborCnt++] = nid;
                            }
                        }
                    }
                }
            }
        }
        neighbors[i][maxNeighborCnt] = neighborCnt;
    }
}

void PBFSolver::solve(Particle* particles) {
    init();
    // apply force & update position
    #pragma omp parallel for
    for (int i = 0; i < para.nParticles; ++i) {
        oldPosition[i] = particles[i].position;
        particles[i].velocity += para.dt * glm::vec3(0.0, -9.8, 0.0);
        particles[i].position += para.dt * particles[i].velocity;
        BoundParticle(particles[i]);
    }
    // find neighboring particles
    buildNeiborhood(particles, para.nParticles, para.h);
    float invH = 1.0 / para.h;
    int times = 0;
    while (times++ < para.iterations)
    {
        // calculate lambda
        #pragma parllel for
        for (int i = 0; i < para.nParticles; ++i) {
            // estimate rho
            float rho = 0;
            float gradMagSum = 0;
            glm::vec3 gradSum(0.0);
            for (int j = 0; j < neighbors[i][maxNeighborCnt]; ++j) {
                int nid = neighbors[i][j];
                rho += Wpoly6(particles[i].position - particles[nid].position, para.h);
                glm::vec3 grad = delWspiky(particles[i].position - particles[nid].position, para.h);
                gradSum += grad;
                gradMagSum += glm::dot(grad, grad);
            }
            float Ci = rho / para.restRho - 1;
            gradMagSum += glm::dot(gradSum, gradSum);
            gradMagSum *= 1.0 / (para.restRho * para.restRho);
            lambda[i] = -Ci / (gradMagSum + para.epsilon);
            printf("rho: %f\n", rho);
        }
        #pragma parallel for 
        for (int i = 0; i < para.nParticles; ++i) {
            glm::vec3 offset(0.0);

            for (int j = 0; j < neighbors[i][maxNeighborCnt]; ++j) {
                int nid = neighbors[i][j];
                float kernelVal = Wpoly6(particles[i].position - particles[nid].position, para.h);
                glm::vec3 kernelGrad = delWspiky(particles[i].position - particles[nid].position, para.h);
                offset += (lambda[i] + lambda[nid] + Scorr(kernelVal)) * kernelGrad;
            }
            offset /= para.restRho;
            particlesOffset[i] = offset;
        }

        #pragma parallel for 
        for (int i = 0; i < para.nParticles; ++i) {
            // collsion detection
            particles[i].position += particlesOffset[i];
            collsionProcess(particles[i]);
        }
    }
    #pragma parallel for 
    for (int i = 0; i < para.nParticles; ++i) {
        particles[i].velocity = (particles[i].position - oldPosition[i]) * para.invDt;
    }
    //#pragma omp parallel for 
    for (int i = 0; i < para.nParticles; ++i) {
        // update velocity
        glm::vec3 curl(0.0);

        for (int j = 0; j < neighbors[i][maxNeighborCnt]; ++j) {
            int nid = neighbors[i][j];
            glm::vec3 vij = particles[i].velocity - particles[nid].velocity;
            curl += glm::cross(vij, delWspiky(particles[i].position - particles[nid].position, para.h));
        }
        particlesCurl[i] = curl;
    }
    #pragma omp parallel for
    for (int i = 0; i < para.nParticles; ++i) {
        glm::vec3 viscosity(0.0);
        glm::vec3 eta(0.0);
        
        for (int j = 0; j < neighbors[i][maxNeighborCnt]; ++j) {
            int nid = neighbors[i][j];
            glm::vec3 vij = particles[i].velocity - particles[nid].velocity;
            // viscosity
            viscosity -= vij * Wpoly6(particles[i].position - particles[nid].position, para.h) ;
            eta += glm::length(particlesCurl[nid]) * Wpoly6(particles[i].position - particles[nid].position, para.h);
        }
        if (!equalZero(glm::length(eta))) {
            eta = glm::normalize(eta);
        }
        // apply vorticity comfinement & XSPH viscosity
        particlesVorticity[i] = fVorticity(eta, particlesCurl[i]);
        particlesViscosity[i] = para.c * viscosity;
    }
    #pragma omp parallel for
    for (int i = 0; i < para.nParticles; ++i) {
        particles[i].velocity += para.dt * particlesVorticity[i] + particlesViscosity[i];
    }

}