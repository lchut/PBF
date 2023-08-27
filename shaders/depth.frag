#version 450 core

layout(location = 0) out float depth;

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

    vec3 surfacePos = viewPos.xyz + radius * vNormal;
    vec4 projSurfacePos = projection * vec4(surfacePos, 1.0f);
    gl_FragDepth =  0.5 * (projSurfacePos.z / projSurfacePos.w) + 0.5;
    depth = -surfacePos.z;
}