#include "AABB.h"
#include <algorithm>

AABB Union(const AABB& box1, const AABB& box2) {
    glm::vec3 pMin = glm::vec3(
        std::min(box1.pMin.x, box2.pMin.x),
        std::min(box1.pMin.y, box2.pMin.y),
        std::min(box1.pMin.z, box2.pMin.z)
    );
    glm::vec3 pMax = glm::vec3(
        std::max(box1.pMin.x, box2.pMin.x),
        std::max(box1.pMin.y, box2.pMin.y),
        std::max(box1.pMin.z, box2.pMin.z)
    );
    return AABB(pMin, pMax);
}

bool inside(const AABB& box, const glm::vec3& p) {
    return !(p.x < box.pMin.x || p.x > box.pMax.x ||
        p.y < box.pMin.y || p.y > box.pMax.y ||
        p.z < box.pMin.z || p.z > box.pMax.z);
}