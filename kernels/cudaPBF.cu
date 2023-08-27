#include "cudaPBF.h"
#include "cudaPBFSolver.cuh"
#include <chrono>
void initCUDASolverDataBuffer(const std::vector<glm::vec3>& position, const std::vector<glm::vec3>& velocity,
    SolverData& data, const SolverParameter& para,
    const AABB& gridSpace, const AABB& boundary) {
    //CHECK_CUDA(cudaMalloc((void**)&data.position, para.nParticles * 3 * sizeof(float)));
    CHECK_CUDA(cudaMalloc((void**)&data.velocity, para.nParticles * 3 * sizeof(float)));
    CHECK_CUDA(cudaMalloc((void**)&data.oldPosition, para.nParticles * 3 * sizeof(float)));
    CHECK_CUDA(cudaMalloc((void**)&data.lambda, para.nParticles * sizeof(float)));
    CHECK_CUDA(cudaMalloc((void**)&data.particlesOffset, para.nParticles * 3 * sizeof(float)));
    CHECK_CUDA(cudaMalloc((void**)&data.particlesCurl, para.nParticles * 3 * sizeof(float)));
    CHECK_CUDA(cudaMalloc((void**)&data.particlesVorticity, para.nParticles * 3 * sizeof(float)));
    CHECK_CUDA(cudaMalloc((void**)&data.particlesViscosity, para.nParticles * 3 * sizeof(float)));
    CHECK_CUDA(cudaMalloc((void**)&data.neighbors, para.nParticles * (MAX_NEIGHBOR + 1) * sizeof(int)));
    CHECK_CUDA(cudaMalloc((void**)&data.gridSpace, 6 * sizeof(float)));
    CHECK_CUDA(cudaMalloc((void**)&data.boundary, 6 * sizeof(float)));
    unsigned int gridCnt = UNIFORM_GRID_DIM_SIZE_X * UNIFORM_GRID_DIM_SIZE_Y * UNIFORM_GRID_DIM_SIZE_Z;
    unsigned int blockCnt = gridCnt / (PREFIX_BLOCK_DIM_X * PREFIX_BLOCK_DIM_Y);
    CHECK_CUDA(cudaMalloc((void**)&data.grid, gridCnt * sizeof(Grid_)));
    CHECK_CUDA(cudaMalloc((void**)&data.prefixSum, gridCnt * sizeof(int)));
    CHECK_CUDA(cudaMalloc((void**)&data.blockSum, blockCnt * sizeof(int)));
    CHECK_CUDA(cudaMalloc((void**)&data.gridParticlesCnt, gridCnt * sizeof(int)));
    CHECK_CUDA(cudaMalloc((void**)&data.gridParticlesID, para.nParticles * sizeof(int)));

    cudaMemcpy(data.position, position.data(), para.nParticles * 3 * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(data.velocity, velocity.data(), para.nParticles * 3 * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(data.gridSpace, &gridSpace, 6 * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(data.boundary, &boundary, 6 * sizeof(float), cudaMemcpyHostToDevice);
}

void freeCUDASolverDataBuffer(SolverData& data) {
    CHECK_CUDA(cudaFree(data.velocity));
    CHECK_CUDA(cudaFree(data.oldPosition));
    CHECK_CUDA(cudaFree(data.lambda));
    CHECK_CUDA(cudaFree(data.particlesOffset));
    CHECK_CUDA(cudaFree(data.particlesCurl));
    CHECK_CUDA(cudaFree(data.particlesVorticity));
    CHECK_CUDA(cudaFree(data.particlesViscosity));
    CHECK_CUDA(cudaFree(data.neighbors));
    CHECK_CUDA(cudaFree(data.gridSpace));
    CHECK_CUDA(cudaFree(data.boundary));
    CHECK_CUDA(cudaFree(data.grid));
    CHECK_CUDA(cudaFree(data.prefixSum));
    CHECK_CUDA(cudaFree(data.blockSum));
    CHECK_CUDA(cudaFree(data.gridParticlesCnt));
    CHECK_CUDA(cudaFree(data.gridParticlesID));
}

void cudaPBFSolve(SolverData& data, const SolverParameter& para) {

    unsigned int unifromGridCnt = UNIFORM_GRID_DIM_SIZE_X * UNIFORM_GRID_DIM_SIZE_Y * UNIFORM_GRID_DIM_SIZE_Z;
    cudaMemset(data.gridParticlesCnt, 0, sizeof(int) * unifromGridCnt);
    int N = para.nParticles;
    dim3 block(BLOCK_DIM_X, BLOCK_DIM_Y);
    int blockSize = BLOCK_DIM_X * BLOCK_DIM_Y;
    dim3 grid((N + blockSize - 1) / blockSize, 1);
    // update position & velocity
    int* temp = new int[unifromGridCnt];
    int* gpuResult = new int[unifromGridCnt];

    updateParticles<<<grid, block>>>(data.position, data.velocity, data.oldPosition, para.nParticles, para.dt, data.boundary);
    CHECK_CUDA(cudaDeviceSynchronize());
    // build neiborhood

    calculateParticlesGid<<<grid, block>>>(data.gridParticlesCnt, data.position, para.nParticles, para.h, data.gridSpace);
    CHECK_CUDA(cudaDeviceSynchronize());

    int prefixBlockSize = PREFIX_BLOCK_DIM_X * PREFIX_BLOCK_DIM_Y;
    int blockCnt = (unifromGridCnt + prefixBlockSize - 1) / prefixBlockSize;
    calculatePrefixSum <<< dim3(blockCnt, 1), dim3(PREFIX_BLOCK_DIM_X >> 1, PREFIX_BLOCK_DIM_Y) >>> (data.gridParticlesCnt, data.prefixSum, data.blockSum);
    CHECK_CUDA(cudaDeviceSynchronize());

    calculatePrefixSum<<<dim3(1, 1), dim3(32, blockCnt / 64)>>>(data.blockSum, data.blockSum, nullptr);
    CHECK_CUDA(cudaDeviceSynchronize());

    addBlockSum<<<dim3(blockCnt, 1), dim3(PREFIX_BLOCK_DIM_X, PREFIX_BLOCK_DIM_Y) >> >(data.prefixSum, data.blockSum, unifromGridCnt);
    CHECK_CUDA(cudaDeviceSynchronize());

    setGrid<<<dim3(blockCnt, 1), dim3(PREFIX_BLOCK_DIM_X, PREFIX_BLOCK_DIM_Y) >>>(data.grid, data.gridParticlesCnt, data.prefixSum, unifromGridCnt);
    CHECK_CUDA(cudaDeviceSynchronize());

    countingSort<<<grid, block>>>(data.gridParticlesID, data.gridParticlesCnt, data.grid, data.position, para.nParticles, para.h, data.gridSpace);
    CHECK_CUDA(cudaDeviceSynchronize());

    buildNeighborhood<<<grid, block>>>(data.neighbors, data.gridParticlesID, data.position, data.grid, para.nParticles, para.h, data.gridSpace);
    CHECK_CUDA(cudaDeviceSynchronize());

    int itrCnt = 0;
    while(itrCnt++ < para.iterations) {
        // calculate Lambda
        calculateLambda<<<grid, block>>>(data.lambda, data.position, data.neighbors, para.nParticles, para.h, para.restRho, para.epsilon);
        CHECK_CUDA(cudaDeviceSynchronize());
        // calculate offset
        calculateOffset<<<grid, block>>>(data.particlesOffset, data.position, data.lambda, data.neighbors, para);
        CHECK_CUDA(cudaDeviceSynchronize());
        // update position & collision process
        updatePosition<<<grid, block>>>(data.position, data.velocity, data.oldPosition, data.particlesOffset, para.nParticles, para.dt, data.boundary);
        CHECK_CUDA(cudaDeviceSynchronize());
    }

    //
    // calculate curl
    calculateCurl<<<grid, block>>>(data.particlesCurl, data.position, data.velocity, data.neighbors, para.nParticles, para.h);
    CHECK_CUDA(cudaDeviceSynchronize());
    // calculate vorticity & XSPH viscosity
    calculateVorticityAndViscosity<<<grid, block>>>(data.particlesVorticity, data.particlesViscosity,
                                data.position, data.velocity, data.particlesCurl, data.neighbors,
                                para.nParticles, para.h, para.c, para.vorticityConfinementEpilon);
    CHECK_CUDA(cudaDeviceSynchronize());
    // apply vorticity confinement and XSPH viscosity
    applyVorticityAndViscosity<<<grid, block>>>(data.velocity, data.particlesViscosity, data.particlesVorticity, para.nParticles, para.dt);
    CHECK_CUDA(cudaDeviceSynchronize());

}