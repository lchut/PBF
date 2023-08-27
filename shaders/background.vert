#version 450 core

layout(location = 0) in vec3 mPosition;
layout(location = 1) in vec3 mNormal;
layout(location = 2) in vec3 mColor;

out vec3 vColor;
out vec3 posWorld;
out vec3 vNormal;

uniform mat4 MVP;

void main() {
    vColor = mColor;
    posWorld = mPosition;
    vNormal = mNormal;
    gl_Position = MVP * vec4(mPosition, 1.0f);
}