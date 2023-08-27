#version 450 core

layout(location=0) out vec4 fragColor;

in vec2 texCoord;

uniform sampler2D tex;

void main() {
    float r = texture(tex, texCoord).x;
    fragColor = vec4(vec3(r/20), 1.0);
    //fragColor = vec4(0.0, 0.0, 10.0 * texture(tex, texCoord).z, 1.0);
}