#include "test.h"
#include "solver.h"

void testSolver() {
    SolverParameter para;
    para.h = 0.1;
    para.nParticles = 2000;
    para.restRho = 9.0 * 315 * 125 * invPi/ (64 * 729 * para.h * para.h * para.h);
    para.epsilon = 0.0001f;
    para.k = 0.1;
    para.deltaq = 0.1 * para.h;
    para.n = 4;
    para.c = 0.01f;
    para.dt = 0.00016;
    para.invDt = 1.0 / para.dt;
    para.iterations = 4;
    para.box = AABB(glm::vec3(-2.0 , 0.1, -2.0), glm::vec3(2.0, 2.1, 2.0));
    PBFSolver solver(para);
    AABB initRegion(glm::vec3(-1.0, 0.5, -1.0), glm::vec3(1.0, 1.0, 1.0));
    ParticleManager particleManager(para.nParticles);
    auto start = std::chrono::system_clock::now();
    auto end = std::chrono::system_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << double(duration.count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den << std::endl;
}