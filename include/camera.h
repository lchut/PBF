#ifndef PBF_CAMERA_H_
#define PBF_CAMERA_H_

#include "config.h"

class PerspectiveCamera {
public:
    PerspectiveCamera(const glm::vec3& _pos,
    const glm::vec3& _target,
    const glm::vec3& _up) :
    pos(_pos), target(_target), up(_up) {
        in = glm::normalize(pos - target);
        right = glm::normalize(glm::cross(up, in));
        up = glm::normalize(glm::cross(in, right));
    }
    inline glm::mat4 getViewMatrix() const {
        return glm::lookAt(pos, pos - in, up);
    }

    glm::vec3 pos;
    glm::vec3 target;
    glm::vec3 in;
    glm::vec3 up;
    glm::vec3 right;
};

#endif