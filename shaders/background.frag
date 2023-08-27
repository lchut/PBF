#version 450 core

out vec4 fragColor;

in vec3 vColor;
in vec3 posWorld;
in vec3 vNormal;

uniform vec3 camPos;

struct PointLight {
    vec3 pos;
    vec3 diffuse;
    vec3 specular;
};

uniform PointLight pointLight;

void main() {
    vec3 lightDir = normalize(pointLight.pos - posWorld);
    vec3 viewDir = normalize(camPos - posWorld);

    vec3 ambient = vec3(0.2, 0.2, 0.2);
    float diff = max(dot(lightDir, vNormal), 0.0);
    vec3 diffuse = diff * pointLight.diffuse * vColor;
    
    fragColor = vec4(ambient + diffuse, 1.0f);
}

