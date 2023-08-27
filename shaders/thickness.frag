#version 450 core

layout(location = 0) out float thickness;

in vec4 viewPos;

uniform mat4 projection;
uniform mat4 view;
uniform float radius;

void main() {
    vec3 vNormal;
    vNormal.xy = vec2(2.0, -2.0) * gl_PointCoord + vec2(-1.0, 1.0);
    float posDis = vNormal.x * vNormal.x + vNormal.y * vNormal.y;
    if (posDis > 1) {
        discard;
        return;
    }
    vNormal.z = sqrt(1 - posDis);

    thickness = 2.0 * vNormal.z * radius;
}