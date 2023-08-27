#version 450 core

layout(location=0) out vec3 normal;
in vec2 texCoord;

uniform float width;
uniform float height;
uniform float Fx;
uniform float Fy;

uniform sampler2D zValTex;

float sampleZVal(vec2 pos) {
    return texture(zValTex, pos).x;
} 

void main() {
    float Cx = 2.0 / (width * Fx);
    float Cy = 2.0 / (height * Fy);
    float dx = 1.0 / width;
    float dy = 1.0 / height;
    vec3 n = vec3(0.0, 0.0, 0.0);
    float z = sampleZVal(texCoord);
    if (z > 0) {
        float zPlusDx = sampleZVal(texCoord + vec2(dx, 0.0));
        float zMinusDx = sampleZVal(texCoord + vec2(-dx, 0.0));
        float dzdx = (zPlusDx == 0 || zMinusDx == 0) ? 0.0f : (zPlusDx - zMinusDx) * 0.5;

        float zPlusDy = sampleZVal(texCoord + vec2(0.0, dy));
        float zMinusDy = sampleZVal(texCoord + vec2(0.0, -dy));
        float dzdy = (zPlusDy == 0 || zMinusDy == 0) ? 0.0f : (zPlusDy - zMinusDy) * 0.5;

        n = vec3(-Cy * dzdx, Cx * dzdy, Cx*Cy*z);
        float nLen = sqrt(dot(n, n));
        n /= nLen;
    }
    normal = n;
}