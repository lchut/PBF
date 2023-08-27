#include <cuda_runtime.h>
#include <math.h>
#include <iostream>
#include "solverBase.h"
inline void CHECK_CUDA(cudaError_t err) {
    if(err != cudaSuccess) {
        std::cerr << "Error: " << cudaGetErrorString(err) << std::endl;
        exit(-1);
    }
}

struct Vec3 {
    float x, y, z;
    __device__ Vec3(float v) : x(v), y(v), z(v) {}
    __device__ Vec3(float _x, float _y, float _z) : x(_x), y(_y), z(_z) {}
    __inline__ __device__ float length() {
        return sqrtf(x*x + y*y + z*z);
    }
    __inline__ __device__ float lengthSquare() {
        return x*x + y*y + z*z;
    }
    __inline__ __device__ Vec3 operator*(float u) const {
        return Vec3(x * u, y * u, z * u);
    }
    __inline__ __device__ Vec3 operator/(float u) const {
        return Vec3(x / u, y / u, z / u);
    }
    __inline__ __device__ Vec3 operator+(const Vec3& v) const {
        return Vec3(x + v.x, y + v.y, z + v.z);
    }
    __inline__ __device__ Vec3 operator-(const Vec3& v) const {
        return Vec3(x - v.x, y - v.y, z - v.z);
    }
    __inline__ __device__ Vec3& operator+=(const Vec3& v) {
        x += v.x;
        y += v.y;
        z += v.z;
        return *this;
    }
    __inline__ __device__ Vec3& operator-=(const Vec3& v) {
        x -= v.x;
        y -= v.y;
        z -= v.z;
        return *this;
    }
    __inline__ __device__ Vec3& operator/=(float u) {
        x /= u;
        y /= u;
        z /= u;
        return *this;
    }
};

struct AABB_ {
    Vec3 pMin, pMax;
    __device__ AABB_() : pMin(Vec3(1.0)), pMax(Vec3(-1.0)) {}
    __device__ AABB_(const Vec3& _pMin, const Vec3& _pMax) :
        pMin(_pMin), pMax(_pMax) {}
};


struct Grid_ {
    unsigned int startIdx, endIdx;
    __device__ Grid_() : startIdx(0), endIdx(0) {}
};

__inline__ __device__ Vec3 operator*(float u, const Vec3& v) {
    return Vec3(u * v.x, u * v.y, u * v.z);
}

__inline__ __device__ float length(const Vec3& v) {
    return sqrtf(v.x * v.x + v.y * v.y + v.z * v.z);
}

__inline__ __device__ float lengthSquare(const Vec3& v) {
    return v.x * v.x + v.y * v.y + v.z * v.z;
}

__inline__ __device__ Vec3 normalize(const Vec3& v) {
    return v / length(v);
}

__inline__ __device__ float dot(const Vec3& v1, const Vec3& v2) {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

__inline__ __device__ Vec3 cross(const Vec3& v1, const Vec3& v2) {
    return Vec3(v1.y * v2.z - v1.z * v2.y,
        v1.z * v2.x - v1.x * v2.z,
        v1.x * v2.y - v1.y * v2.x);
}

__inline__ __device__ float clamp(float L, float R, float v) {
    if (v < L) { return L; }
    else if (v > R) { return R; }
    return v;
}

__inline__ __device__ float Wpoly6(const Vec3& r, float h) {
    if (dot(r, r) > h * h) { return 0; }
    float c = 315.0f * INVPI / (64.0f * powf(h, 9));
    return c * powf(h * h - dot(r, r), 3);
}

__inline__ __device__ Vec3 WspikyGrad(const Vec3& r, float h) {
    float rLen = length(r);
    if (rLen > h || rLen < 0.01f) { return Vec3(0.0f); }
    float c = -45.0f * INVPI / powf(h, 6);
    return c * powf(h - length(r), 2) * normalize(r);
}

__inline__ __device__ float frac(float v) {
    return v - (int)v;
}

__inline__ __device__ float randomOffset(float u, float v) {
    float noise = frac(sinf(2.0f * (u * 12.9898f + v * 78.233f)) * 43758.5453f);
    return fabsf(noise);
}

__inline__ __device__ void BoundParticle(Vec3& position, const AABB_& boundary) {
    if (position.x < boundary.pMin.x) {
        position.x = boundary.pMin.x + RANDOM_OFFSET_COEFF * randomOffset(position.y, position.z);
    }
    if (position.x > boundary.pMax.x) {
        position.x = boundary.pMax.x - RANDOM_OFFSET_COEFF * randomOffset(position.y, position.z);
    }
    if (position.y < boundary.pMin.y) {
        position.y = boundary.pMin.y + RANDOM_OFFSET_COEFF * randomOffset(position.x, position.z);
    }
    if (position.y > boundary.pMax.y) {
        position.y = boundary.pMax.y - RANDOM_OFFSET_COEFF * randomOffset(position.y, position.z);
    }
    if (position.z < boundary.pMin.z) {
        position.z = boundary.pMin.z + RANDOM_OFFSET_COEFF * randomOffset(position.x, position.y);
    }
    if (position.z > boundary.pMax.z) {
        position.z = boundary.pMax.z - RANDOM_OFFSET_COEFF * randomOffset(position.x, position.y);
    }
}

__inline__ __device__ float Scorr(const Vec3& r, float k, float h, float n, float deltaq) {
    return -k * powf(Wpoly6(r, h) / Wpoly6(Vec3(deltaq, 0.0f, 0.0f), h), n);
}

__inline__ __device__ Vec3 fVorticity(const Vec3& N, const Vec3& w, float epsilon) {
    return epsilon * cross(N, w);
}

__inline__ __device__ bool outOfRange(const int dim[3]) {
    return dim[0] < 0 || dim[0] >= UNIFORM_GRID_DIM_SIZE_X ||
        dim[1] < 0 || dim[1] >= UNIFORM_GRID_DIM_SIZE_Y ||
        dim[2] < 0 || dim[2] >= UNIFORM_GRID_DIM_SIZE_Z;
}

__inline__ __device__ void computeDimID(const Vec3& position, const Vec3& pMin, float invH, int* dimX, int* dimY, int* dimZ) {
    *dimX = int((position.x - pMin.x) * invH);
    *dimY = int((position.y - pMin.y) * invH);
    *dimZ = int((position.z - pMin.z) * invH);
}


__inline__ __device__ int computeGid(const int dim[3]) {
    if (outOfRange(dim)) {
        return -1;
    }
    else {
        return dim[0] + dim[1] * UNIFORM_GRID_DIM_SIZE_X + dim[2] * UNIFORM_GRID_DIM_SIZE_X * UNIFORM_GRID_DIM_SIZE_Y;
    }
}

__inline__ __device__ int computeGid(const Vec3& position, const Vec3& pMin, float invH) {
    int dim[3];
    computeDimID(position, pMin, invH, &dim[0], &dim[1], &dim[2]);
    if (outOfRange(dim)) {
        return -1;
    }
    else {
        return dim[0] + dim[1] * UNIFORM_GRID_DIM_SIZE_X + dim[2] * UNIFORM_GRID_DIM_SIZE_X * UNIFORM_GRID_DIM_SIZE_Y;
    }
}

extern __global__ void updateParticles(Vec3* position, Vec3* velocity, Vec3* oldPosition, int nParticles, float dt, const AABB_* boundary);
extern __global__ void updatePosition(Vec3* position, Vec3* velocity, const Vec3* oldPosition, const Vec3* offset, int nParticles, float dt, const AABB_* boundary);
extern __global__ void applyVorticityAndViscosity(Vec3* velocity, const Vec3* viscosity, const Vec3* vorticity, int nParticles, float dt);
extern __global__ void calculateParticlesGid(int* gridParticlesCnt, const Vec3* position, int nParticles, float h, const AABB_* gridSpace);

extern __global__ void calculatePrefixSum(int* indata, int* outdata, int* blockSum);
extern __global__ void addBlockSum(int* prefixSum, int* blockSum, int N);

extern __global__ void setGrid(Grid_* grid, const int* gridParticlesCnt, const int* prefixSum, int N);
extern __global__ void countingSort(int* gridParticlesID, int* gridParticlesCnt, const Grid_* grid, const Vec3* position, int nParticles, float h, const AABB_* gridSpace);
extern __global__ void buildNeighborhood(int* neighbors, const int* gridParticlesID, const Vec3* position, const Grid_* grid, int nParticles, float h, const AABB_* gridSpace);


extern __global__ void calculateLambda(float* lambda, const Vec3* position, const int* neighbors, int nParticles, float h, float restRho, float epsilon);
extern __global__ void calculateOffset(Vec3* offset, const Vec3* position, const float* lambda, const int* neighbors,
                                const SolverParameter para);
extern __global__ void calculateCurl(Vec3* curl, const Vec3* position, const Vec3* velocity, const int* neighbors,
                                int nParticles, float h);
extern __global__ void calculateVorticityAndViscosity(Vec3* vorticity ,Vec3* viscosity,
                                                const Vec3* position, const Vec3* velocity, const Vec3* curl,
                                                const int* neighbors, int nParticles, float h, float c, float vorticityConfinementEpilon);



