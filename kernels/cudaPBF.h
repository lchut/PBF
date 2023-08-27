#include "solverBase.h"
#include "AABB.h"
#include "particle.h"
#include "cuda_runtime.h"
#include "cuda_gl_interop.h"


void initCUDASolverDataBuffer(const std::vector<glm::vec3>& position, const std::vector<glm::vec3>& velocity,
    SolverData& data, const SolverParameter& para,
    const AABB& gridSpace, const AABB& boundary);
void freeCUDASolverDataBuffer(SolverData& data);
void cudaPBFSolve(SolverData& data, const SolverParameter& para);