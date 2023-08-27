#version 450 core

out vec4 fragColor;

in vec4 viewPos;

uniform mat4 projection;
uniform mat4 view;
uniform float radius;

uniform vec3 camPos;

struct PointLight {
    vec3 pos;
    vec3 diffuse;
    vec3 specular;
};

uniform PointLight pointLight;

void main() {
    vec3 vNormal;
    vNormal.xy = vec2(2.0, -2.0) * gl_PointCoord + vec2(-1.0, 1.0);
    float posDis = vNormal.x * vNormal.x + vNormal.y * vNormal.y;
    if (posDis > 1) discard;
    vNormal.z = sqrt(1 - posDis);
    vec3 lightPos = (view * vec4(pointLight.pos, 1.0f)).xyz;
    vec3 surfacePos = viewPos.xyz + radius * vNormal;
    vec3 lightDir = normalize(lightPos - surfacePos);
    vec3 ambient = vec3(0.2, 0.2, 0.2);
    float diff = max(dot(lightDir, vNormal), 0.0);
    vec3 diffuse = diff * pointLight.diffuse * vec3(0.156, 0.376, 0.705);
    fragColor = vec4(ambient + diffuse, 1.0f);
}