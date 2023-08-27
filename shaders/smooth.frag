#version 450 core

layout(location = 0) out float smoothDepth;

in vec2 texCoord;

uniform float width;
uniform float height;
uniform float Fx;
uniform float Fy;
uniform float smoothingDt;

uniform sampler2D zValTex;

float sampleZVal(vec2 pos) {
    return texture(zValTex, pos).x;
}

void main() {
    float Cx = 2.0f / (width * Fx);
    float Cy = 2.0f / (height * Fy);
    float dx = 1.0f / width;
    float dy = 1.0f / height;
    
    float z = sampleZVal(texCoord);
    if (z > 0) {
        float zPlusDx = sampleZVal(texCoord + vec2(dx, 0.0));
        float zMinusDx = sampleZVal(texCoord + vec2(-dx, 0.0));
        float dzdx = (zPlusDx == 0 || zMinusDx == 0) ? 0.0f : (zPlusDx - zMinusDx) * 0.5;
        float dz2dx2 = (zPlusDx + zMinusDx - 2.0f * z);

        float zPlusDy = sampleZVal(texCoord + vec2(0.0, dy));
        float zMinusDy = sampleZVal(texCoord + vec2(0.0, -dy));
        float dzdy = (zPlusDy == 0 || zMinusDy == 0) ? 0.0f : (zPlusDy - zMinusDy) * 0.5;
        float dz2dy2 = (zPlusDy + zMinusDy - 2.0f * z);

        float zPlusDxPlusDy = sampleZVal(texCoord + vec2(dx, dy));
        float zPlusDxMinusDy = sampleZVal(texCoord + vec2(dx, -dy));
        float zMinusDxPlusDy = sampleZVal(texCoord + vec2(-dx, dy));
        float zMinusDxMinusDy = sampleZVal(texCoord + vec2(-dx, -dy));
        float dz2dxdy = (zPlusDxPlusDy + zMinusDxMinusDy - zMinusDxPlusDy - zPlusDxMinusDy) * 0.25;
        
        vec3 n = vec3(-Cy * dzdx, Cx * dzdy, Cx*Cy*z);
        float D = dot(n, n);
        float Cx2 = Cx * Cx;
        float Cy2 = Cy * Cy;
        float dDdx = 2.0f * (Cy2 * dzdx * dz2dx2 + Cx2 * dzdy * dz2dxdy + Cx2 * Cy2 * z * dzdx);
        float dDdy = 2.0f * (Cy2 * dzdx * dz2dxdy + Cx2 * dzdy * dz2dy2 + Cx2 * Cy2 * z * dzdy);
        float Ex = 0.5f * dzdx * dDdx - dz2dx2 * D;
        float Ey = 0.5f * dzdy * dDdy - dz2dy2 * D;
        float H = 0.5f * (Cy * Ex + Cx * Ey) / (D * sqrt(D));

        smoothDepth = z - H * 0.01;

        if (smoothDepth < 0)
        {
            smoothDepth = z;
        }
    }
    else {
        smoothDepth = 0.0;
    }
    
}