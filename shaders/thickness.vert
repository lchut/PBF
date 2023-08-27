#version 450 core

layout(location = 0) in vec3 mPosition;

out vec4 viewPos;

uniform mat4 projection;
uniform mat4 view;
uniform float radius;
uniform float pointCoeff;

void main() {
    viewPos = view * vec4(mPosition, 1.0f);
    gl_Position = projection * viewPos;
    gl_PointSize = -pointCoeff * radius / viewPos.z; 
}