#include "cudaPBFSolver.cuh"
#include <stdio.h>

__global__ void updateParticles(Vec3* position, Vec3* velocity, Vec3* oldPosition, int nParticles, float dt, const AABB_* boundary) {
    unsigned int blockId = blockIdx.y * gridDim.x + blockIdx.x;
    unsigned int idx = blockId * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x + threadIdx.x);
    if (idx >= nParticles) { return; }
    oldPosition[idx] = position[idx];
    velocity[idx] += dt * Vec3(0.0f, -49.0f, 0.0f);
    position[idx] += dt * velocity[idx];
    BoundParticle(position[idx], *boundary);

}

__global__ void updatePosition(Vec3* position, Vec3* velocity, const Vec3* oldPosition, const Vec3* offset, int nParticles, float dt, const AABB_* boundary) {
    unsigned int blockId = blockIdx.y * gridDim.x + blockIdx.x;
    unsigned int idx = blockId * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x + threadIdx.x);
    if (idx >= nParticles) { return; }
    position[idx] += offset[idx];
    BoundParticle(position[idx], *boundary);
    velocity[idx] = (position[idx] - oldPosition[idx]) / dt;
}

__global__ void applyVorticityAndViscosity(Vec3* velocity, const Vec3* viscosity, const Vec3* vorticity, int nParticles, float dt) {
    unsigned int blockId = blockIdx.y * gridDim.x + blockIdx.x;
    unsigned int idx = blockId * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x + threadIdx.x);
    if (idx >= nParticles) { return; }
    velocity[idx] += dt * vorticity[idx] + viscosity[idx];
}

__global__ void calculateParticlesGid(int* gridParticlesCnt, const Vec3* position, int nParticles, float h, const AABB_* gridSpace) {
    unsigned int blockId = blockIdx.y * gridDim.x + blockIdx.x;
    unsigned int idx = blockId * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x + threadIdx.x);
    if (idx >= nParticles) { return; }
    float invH = 1.0f / h;
    int gid = computeGid(position[idx], gridSpace->pMin, invH);
    if (gid >= 0) {
        atomicAdd(&gridParticlesCnt[gid], 1);
    }
}

__global__ void calculatePrefixSum(int* indata, int* outdata, int* blockSum) {
    __shared__ int temp[SHARED_MEMORY_SIZE];
    unsigned int blockId = blockIdx.y * gridDim.x + blockIdx.x;
    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;
    unsigned int idx = blockId * (blockDim.x * blockDim.y) + tid;
    if (tid == 0) {
        for (int i = 0; i < SHARED_MEMORY_SIZE; ++i) {
            temp[i] = 0;
        }
    }
    __syncthreads();
    int n = blockDim.x * blockDim.y * 2;
    temp[2 * tid] = indata[2 * idx];
    temp[2 * tid + 1] = indata[2 * idx + 1];
    int offset = 1;
    for (int d = n >> 1; d > 0; d >>= 1) {
        __syncthreads();
        if (tid < d) {
            int ai = offset * (2 * tid + 1) - 1;
            int bi = offset * (2 * tid + 2) - 1;
            temp[bi] += temp[ai];
        }
        offset *= 2;
    }
    __syncthreads();
    int blockTotalSum = 0;
    if (tid == 0) {
        blockTotalSum = temp[n - 1];
        temp[n - 1] = 0;
    }
    for (int d = 1; d < n; d *= 2) {
        __syncthreads();
        offset >>= 1;
        if (tid < d) {
            int ai = offset * (2 * tid + 1) - 1;
            int bi = offset * (2 * tid + 2) - 1;
            int t = temp[ai];
            temp[ai] = temp[bi];
            temp[bi] += t;
        }
    }
    __syncthreads();
    outdata[2 * idx] = temp[2 * tid];
    outdata[2 * idx + 1] = temp[2 * tid + 1];
    if (tid == 0 && blockSum != nullptr) {
        blockSum[blockId] = blockTotalSum;
    }
}

__global__ void addBlockSum(int* prefixSum, int* blockSum, int N) {
    unsigned int blockId = blockIdx.y * gridDim.x + blockIdx.x;
    unsigned int idx = blockId * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x + threadIdx.x);
    if (idx >= N) { return; }
    prefixSum[idx] += blockSum[blockId];
}


__global__ void setGrid(Grid_* grid, const int* gridParticlesCnt, const int* prefixSum, int N) {
    unsigned int blockId = blockIdx.y * gridDim.x + blockIdx.x;
    unsigned int idx = blockId * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x + threadIdx.x);
    if (idx >= N) { return; }
    if (idx == N-1) {
        grid[idx].startIdx = prefixSum[idx];
        grid[idx].endIdx = prefixSum[idx] + gridParticlesCnt[N-1];
    }
    else {
        grid[idx].startIdx = prefixSum[idx];
        grid[idx].endIdx = prefixSum[idx + 1];
    }
}
__global__ void countingSort(int* gridParticlesID, int* gridParticlesCnt, const Grid_* grid, const Vec3* position, int nParticles, float h, const AABB_* gridSpace) {
    unsigned int blockId = blockIdx.y * gridDim.x + blockIdx.x;
    unsigned int idx = blockId * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x + threadIdx.x);
    if (idx >= nParticles) { return; }
    float invH = 1.0f / h;
    int gid = computeGid(position[idx], gridSpace->pMin, invH);
    if (gid >= 0) {
        int offset = atomicSub(&gridParticlesCnt[gid], 1);
        gridParticlesID[grid[gid].startIdx + offset - 1] = idx;
    }
}

__global__ void buildNeighborhood(int* neighbors, const int* gridParticlesID, const Vec3* position, const Grid_* grid, int nParticles, float h, const AABB_* gridSpace) {
    unsigned int blockId = blockIdx.y * gridDim.x + blockIdx.x;
    unsigned int idx = blockId * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x + threadIdx.x);
    if (idx >= nParticles) { return; }
    int* pNeighbor = &neighbors[idx * (MAX_NEIGHBOR + 1)];
    int neighborCnt = 0;
    int dim[3];
    computeDimID(position[idx], gridSpace->pMin, 1.0f / h, &dim[0], &dim[1], &dim[2]);
    int neiborDim[3];
    for (int j = 0; j < 27 && neighborCnt < MAX_NEIGHBOR; ++j) {
        int xOffset = j % 3 - 1;
        int yOffSet = int(j / 3) % 3 - 1;
        int zOffset = int(j / 9) - 1;
        neiborDim[0] = dim[0] + xOffset;
        neiborDim[1] = dim[1] + yOffSet;
        neiborDim[2] = dim[2] + zOffset;
        int gid = computeGid(neiborDim);
        if (gid >=0) {
            for (int k = grid[gid].startIdx; k < grid[gid].endIdx && neighborCnt < MAX_NEIGHBOR; ++k) {
                int nid = gridParticlesID[k];
                if (nid != idx) {
                    Vec3 r = position[idx] - position[nid];
                    if (dot(r, r) <= h * h) {
                        pNeighbor[neighborCnt++] = nid;
                    }
                }
            }
        }
    }
    pNeighbor[MAX_NEIGHBOR] = neighborCnt;
}

__global__ void calculateLambda(float* lambda, const Vec3* position, const int* neighbors, int nParticles, float h, float restRho, float epsilon) {
    unsigned int blockId = blockIdx.y * gridDim.x + blockIdx.x;
    unsigned int idx = blockId * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x + threadIdx.x);
    if (idx >= nParticles) { return; }
    float rho = 0.0f;
    float gradMagSum = 0;
    Vec3 gradSum(0.0f);
    const int* pNeighbors = &neighbors[idx*(MAX_NEIGHBOR + 1)];
    int neighborsCnt = pNeighbors[MAX_NEIGHBOR];
    for (int j = 0; j < neighborsCnt; ++j) {
        int nid = pNeighbors[j];
        if (nid < 0 || nid >= nParticles) { printf("%d\n", nid); }
        Vec3 r = position[idx] - position[nid];

        rho += Wpoly6(r , h);
        Vec3 grad = WspikyGrad(r, h);
        gradSum += grad;
        gradMagSum += dot(grad, grad);
    }
    //printf("%f\n", rho);
    float Ci = rho / restRho - 1;
    gradMagSum += dot(gradSum, gradSum);
    gradMagSum /= restRho * restRho;
    lambda[idx] = -Ci / (gradMagSum + epsilon);
}

__global__ void calculateOffset(Vec3* offset, const Vec3* position, const float* lambda, const int* neighbors,
                                const SolverParameter para) {
    unsigned int blockId = blockIdx.y * gridDim.x + blockIdx.x;
    unsigned int idx = blockId * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x + threadIdx.x);
    if (idx >= para.nParticles) { return; }
    Vec3 pOffset(0.0f);
    const int* pNeighbors = &neighbors[idx*(MAX_NEIGHBOR + 1)];
    int neighborsCnt = pNeighbors[MAX_NEIGHBOR];
    for (int j = 0; j < neighborsCnt; ++j) {
        int nid = pNeighbors[j];
        Vec3 r = position[idx] - position[nid];
        float kernelVal = Wpoly6(r, para.h);
        Vec3 kernelGrad = WspikyGrad(r, para.h);
        pOffset += (lambda[idx] + lambda[nid] + Scorr(r, para.k, para.h, para.n, para.deltaq)) * kernelGrad;
    }
    pOffset /= para.restRho;
    offset[idx] = pOffset;
}

__global__ void calculateCurl(Vec3* curl, const Vec3* position, const Vec3* velocity, const int* neighbors,
                                int nParticles, float h) {
    unsigned int blockId = blockIdx.y * gridDim.x + blockIdx.x;
    unsigned int idx = blockId * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x + threadIdx.x);
    if (idx >= nParticles) { return; }
    Vec3 pCurl(0.0f);
    const int* pNeighbors = &neighbors[idx*(MAX_NEIGHBOR + 1)];
    int neighborsCnt = pNeighbors[MAX_NEIGHBOR];
    for (int j = 0; j < neighborsCnt; ++j) {
        int nid = pNeighbors[j];
        Vec3 vij = velocity[nid] - velocity[idx];
        pCurl += cross(vij, WspikyGrad(position[idx] - position[nid], h));
    }
    curl[idx] = pCurl;
}

__global__ void calculateVorticityAndViscosity(Vec3* vorticity ,Vec3* viscosity,
                                                const Vec3* position, const Vec3* velocity, const Vec3* curl,
                                                const int* neighbors, int nParticles, float h, float c, float vorticityConfinementEpilon) {
    unsigned int blockId = blockIdx.y * gridDim.x + blockIdx.x;
    unsigned int idx = blockId * (blockDim.x * blockDim.y) + (threadIdx.y * blockDim.x + threadIdx.x);
    if (idx >= nParticles) { return; }
    Vec3 pViscosity(0.0f);
    Vec3 eta(0.0f);
    const int* pNeighbors = &neighbors[idx*(MAX_NEIGHBOR + 1)];
    int neighborsCnt = pNeighbors[MAX_NEIGHBOR];
    for (int j = 0; j < neighborsCnt; ++j) {
        int nid = pNeighbors[j];
        Vec3 vij = velocity[idx] - velocity[nid];
        Vec3 r = position[idx] - position[nid];
        pViscosity -= vij * Wpoly6(r, h);
        eta += length(curl[nid]) * WspikyGrad(position[idx] - position[nid], h);
    }
    if (fabsf(length(eta)) > 1e-3) {
        eta = normalize(eta);
    }
    vorticity[idx] = fVorticity(eta, curl[idx], vorticityConfinementEpilon);
    viscosity[idx] = c * pViscosity;
}
