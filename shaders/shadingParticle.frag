#version 450 core

out vec4 fragColor;
in vec2 texCoord;

struct PointLight {
    vec3 pos;
    vec3 diffuse;
    vec3 specular;
};

uniform PointLight pointLight;

uniform float roi;
uniform float width;
uniform float height;
uniform float tanHalfFov;
uniform mat4 viewMat;

uniform sampler2D zValTex;
uniform sampler2D thicknessTex;
uniform sampler2D normalTex;
uniform sampler2D backgroundTex;

vec3 getSurfacePos() {
    float z = texture(zValTex, texCoord).x;
    float x = (2 * texCoord.x - 1) * tanHalfFov * z;
    float y = (2 * texCoord.y - 1) * tanHalfFov * z;
    return vec3(x, y, -z);
}

float schlick(float cosine) {
    float f = (1 - roi) / (1 + roi);
    float F = f * f;
    return F * (1 - F) * pow(1.0 - cosine, 5.0);
}

void main() {
    vec3 n = texture(normalTex, texCoord).xyz;
    vec3 surfacePos = getSurfacePos();
    vec3 v = normalize(-surfacePos);
    vec3 lightPos = (viewMat * vec4(pointLight.pos, 1.0f)).xyz;
    vec3 lightDir = normalize(lightPos- surfacePos);
    vec3 h = normalize(n + v);

    vec3 Cfliud = vec3(0.1, 0.4, 0.8);
    float thickness = 0.08 * texture(thicknessTex, texCoord).x;
    float beta = 0.01 * thickness;
    vec3 Sxy = texture(backgroundTex, texCoord + beta * vec2(n.x, n.y)).rgb;
    vec3 a = mix(Cfliud, Sxy, exp(-thickness));
    
    float F = schlick(clamp(dot(n, v), 0.0, 1.0));
    
    float diff = max(dot(lightDir, n), 0.0);
    vec3 diffuse = diff * pointLight.diffuse * Cfliud;
    vec3 reflectDir = reflect(-lightDir, n);
    float spec = pow(max(dot(v, reflectDir), 0.0), 32);
    vec3 specular = spec * vec3(1.0, 1.0, 1.0);

    fragColor = vec4(a * (1 - F) + diffuse * F + specular, 1.0);

}