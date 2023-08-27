#ifndef PBF_AABB_H_
#define PBF_AABB_H_

#include "config.h"

class AABB {
public:
    AABB() : pMin(glm::vec3(1.0, 1.0, 1.0)), pMax(glm::vec3(-1.0, -1.0, -1.0)) {}
    AABB(const glm::vec3& p) : pMin(p), pMax(p) {}
    AABB(const glm::vec3& _pMin, const glm::vec3& _pMax) : pMin(_pMin), pMax(_pMax) {}
    glm::vec3 pMin;
    glm::vec3 pMax;
};
AABB Union(const AABB& box1, const AABB& box2);
bool inside(const AABB& box, const glm::vec3& p);
#endif