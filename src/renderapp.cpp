#include "renderapp.h"
#include "camera.h"
#include "program.h"
#include "solver.h"

#include "cudaPBF.h"
#include <ctime>
#include <chrono>

int fluidSceneIdx = 1;
int fluidRenderMode = 0;

void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    glViewport(0, 0, width, height);
}

void RenderApp::run() {
    initWindow();
    initOpenGL();
    mainLoop();
    cleanUp();
}

void RenderApp::initWindow() {
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 4);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    window = glfwCreateWindow(WINDOW_WIDTH, WINDOW_HEIGHT, "PBF", NULL, NULL);
    if (window == NULL)
    {
        glfwTerminate();
        throw std::runtime_error("Failed to create GLFW window");
    }
    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
}

void RenderApp::initOpenGL() {
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        throw std::runtime_error("Failed to initialize GLAD");
    }
}

void RenderApp::mainLoop() {
    const float fov = 45.0f;
    const float nearPlane = 0.1f;
    const float farPlane = 1000.0f;
    const float aspectRatio = float(WINDOW_WIDTH) / WINDOW_HEIGHT;
    PerspectiveCamera camera(glm::vec3(0.0f, 40.0f, 85.0f), glm::vec3(0.0f, 20.0f, -1.0f), glm::vec3(0.0f, 1.0, 0.0));
    glm::mat4x4 model = glm::mat4(1.0f);
    glm::mat4x4 projection = glm::perspective(fov, aspectRatio, nearPlane, farPlane);
    glm::mat4x4 view = camera.getViewMatrix();

    Program backgroundProgram("../shaders/background.vert", "../shaders/background.frag");
    float background[] = {
        // left
        -30.0f, 50.0f, 30.0f,   1.0f, 0.0f, 0.0f, 0.95f, 0.95f, 0.95f,
        -30.0f, 0.0f, 30.0f,  1.0f, 0.0f, 0.0f, 0.95f, 0.95f, 0.95f,
        -30.0f, 0.0f, -30.0f, 1.0f, 0.0f, 0.0f, 0.95f, 0.95f, 0.95f,

        -30.0f, 50.0f, 30.0f,   1.0f, 0.0f, 0.0f, 0.95f, 0.95f, 0.95f,
        -30.0f, 0.0f, -30.0f, 1.0f, 0.0f, 0.0f, 0.95f, 0.95f, 0.95f,
        -30.0f, 50.0f, -30.0f,  1.0f, 0.0f, 0.0f, 0.95f, 0.95f, 0.95f,
        // right
        30.0f, 50.0f, 30.0f,   -1.0f, 0.0f, 0.0f, 0.95f, 0.95f, 0.95f,
        30.0f, 0.0f, -30.0f, -1.0f, 0.0f, 0.0f, 0.95f, 0.95f, 0.95f,
        30.0f, 0.0f, 30.0f,  -1.0f, 0.0f, 0.0f, 0.95f, 0.95f, 0.95f,

        30.0f, 50.0f, 30.0f,   -1.0f, 0.0f, 0.0f, 0.95f, 0.95f, 0.95f,
        30.0f, 50.0f, -30.0f,  -1.0f, 0.0f, 0.0f, 0.95f, 0.95f, 0.95f,
        30.0f, 0.0f, -30.0f, -1.0f, 0.0f, 0.0f, 0.95f, 0.95f, 0.95f,
        // bottom
        -30.0f, 0.0f, -30.0f,  0.0f, 1.0f, 0.0f, 0.95f, 0.95f, 0.95f,
        30.0f, 0.0f, -30.0f,   0.0f, 1.0f, 0.0f, 0.95f, 0.95f, 0.95f,
        30.0f, 0.0f, 30.0f,  0.0f, 1.0f, 0.0f, 0.95f, 0.95f, 0.95f,

        -30.0f, 0.0f, 30.0f,  0.0f, 1.0f, 0.0f, 0.95f, 0.95f, 0.95f,
        30.0f, 0.0f, 30.0f,  0.0f, 1.0f, 0.0f, 0.95f, 0.95f, 0.95f,
        -30.0f, 0.0f, -30.0f, 0.0f, 1.0f, 0.0f, 0.95f, 0.95f, 0.95f,
        // back
        -30.0f, 0.0f, -30.0f, 0.0f, 0.0f, 1.0f, 0.95f, 0.95f, 0.95f,
        30.0f, 0.0f, -30.0f,  0.0f, 0.0f, 1.0f, 0.95f, 0.95f, 0.95f,
        30.0f, 50.0f, -30.0f,   0.0f, 0.0f, 1.0f, 0.95f, 0.95f, 0.95f,

        -30.0f, 0.0f, -30.0f, 0.0f, 0.0f, 1.0f, 0.95f, 0.95f, 0.95f,
        30.0f, 50.0f, -30.0f,   0.0f, 0.0f, 1.0f, 0.95f, 0.95f, 0.9f,
        -30.0f, 50.0f, -30.0f,  0.0f, 0.0f, 1.0f, 0.95f, 0.95f, 0.9f
    };
    uint32_t backgroundVBO;
    uint32_t backgroundVAO;
    glGenVertexArrays(1, &backgroundVAO);
    glGenBuffers(1, &backgroundVBO);
    glBindBuffer(GL_ARRAY_BUFFER, backgroundVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(background), background, GL_STATIC_DRAW);

    glBindVertexArray(backgroundVAO);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 9 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 9 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 9 * sizeof(float), (void*)(6 * sizeof(float)));
    glEnableVertexAttribArray(2);

    backgroundProgram.use();
    backgroundProgram.setMat4("MVP", projection * camera.getViewMatrix() * model);
    backgroundProgram.setVec3("camPos", camera.pos);
    backgroundProgram.setVec3("pointLight.pos", glm::vec3(0.0, 50.0, 00.0));
    backgroundProgram.setVec3("pointLight.diffuse", glm::vec3(0.8, 0.8, 0.8));
    backgroundProgram.setVec3("pointLight.specular", glm::vec3(1.0, 1.0, 1.0));
    Program particleProgram("../shaders/particle.vert", "../shaders/particle.frag");

    int dimParticlesCnt[3] = { UNIFORM_GRID_DIM_SIZE_X, UNIFORM_GRID_DIM_SIZE_Y, UNIFORM_GRID_DIM_SIZE_Z};
    SolverParameter para;
    para.h = 1.0;
    para.nParticles = dimParticlesCnt[0] * dimParticlesCnt[1] * dimParticlesCnt[2];
    para.restRho = 6;
    para.epsilon = 1.0f;
    para.k = 0.01;
    para.deltaq = 0.1 * para.h;
    para.n = 4;
    para.c = 0.1f;
    para.dt = 0.016f;
    para.invDt = 1.0f / para.dt;
    para.iterations = 4;
    para.vorticityConfinementEpilon = 0.1;
    para.radius = 0.3;
    AABB box = AABB(glm::vec3(-30.0 , 0.0, -30.0), glm::vec3(30.0, 40.0, 30.0));
    AABB boundary = AABB(box.pMin + glm::vec3(para.radius), box.pMax - glm::vec3(para.radius));
    AABB initRegion;
    float dimGap[3];
    std::vector<glm::vec3> position(para.nParticles);
    std::vector<glm::vec3> velocity(para.nParticles);
    srand((unsigned int)time(NULL));
    if (fluidSceneIdx == 1) {
        initRegion = AABB(glm::vec3(-16.0, 0.2, -16.0), glm::vec3(16.0, 0.2 + 32, 16.0));
        dimGap[0] = (initRegion.pMax.x - initRegion.pMin.x) / dimParticlesCnt[0];
        dimGap[1] = (initRegion.pMax.y - initRegion.pMin.y) / dimParticlesCnt[1];
        dimGap[2] = (initRegion.pMax.z - initRegion.pMin.z) / dimParticlesCnt[2];
        for (int i = 0; i < para.nParticles; ++i) {
            int dim[3] = {
                i % dimParticlesCnt[0],
                i / dimParticlesCnt[0] % dimParticlesCnt[1],
                i / (dimParticlesCnt[0] * dimParticlesCnt[1])
            };
            position[i] = initRegion.pMin + glm::vec3(
                dim[0] * dimGap[0], dim[1] * dimGap[1], dim[2] * dimGap[2]
            ) + 0.2f * glm::vec3(float(rand())/RAND_MAX, float(rand())/RAND_MAX, float(rand())/RAND_MAX);
            velocity[i] = glm::vec3(0.0, 20.0, 0.0);
        }
    }
    else if (fluidSceneIdx == 2) {
        initRegion = AABB(glm::vec3(-28.0, 0.3, -4.0), glm::vec3(4.0, 0.3 + 32, 28.0));
        dimGap[0] = (initRegion.pMax.x - initRegion.pMin.x) / dimParticlesCnt[0];
        dimGap[1] = (initRegion.pMax.y - initRegion.pMin.y) / dimParticlesCnt[1];
        dimGap[2] = (initRegion.pMax.z - initRegion.pMin.z) / dimParticlesCnt[2];
        for (int i = 0; i < para.nParticles; ++i) {
            int dim[3] = {
                i % dimParticlesCnt[0],
                i / dimParticlesCnt[0] % dimParticlesCnt[1],
                i / (dimParticlesCnt[0] * dimParticlesCnt[1])
            };
            position[i] = initRegion.pMin + glm::vec3(
                dim[0] * dimGap[0], dim[1] * dimGap[1], dim[2] * dimGap[2]
            ) + 0.2f * glm::vec3(float(rand())/RAND_MAX, float(rand())/RAND_MAX, float(rand())/RAND_MAX);
            velocity[i] = glm::vec3(0.0, 0.0, -20.0);
        }
    }
    else if (fluidSceneIdx == 3) {
        int halfN = para.nParticles / 2;
        initRegion = AABB(glm::vec3(-29.7, 0.2, -29.7), glm::vec3(-2.0, 0.2 + 32, -2.0));
        dimGap[0] = (initRegion.pMax.x - initRegion.pMin.x) / dimParticlesCnt[0];
        dimGap[1] = (initRegion.pMax.y - initRegion.pMin.y) / dimParticlesCnt[1];
        dimGap[2] = (initRegion.pMax.z - initRegion.pMin.z) / dimParticlesCnt[2];
        for (int i = 0; i < halfN; ++i) {
            int dim[3] = {
                i % dimParticlesCnt[0],
                i / dimParticlesCnt[0] % 32,
                i / (dimParticlesCnt[0] * 32)
            };
            position[i] = initRegion.pMin + glm::vec3(
                dim[0] * dimGap[0], dim[1] * dimGap[1], dim[2] * dimGap[2]
            ) + 0.2f * glm::vec3(float(rand())/RAND_MAX, float(rand())/RAND_MAX, float(rand())/RAND_MAX);
            velocity[i] = glm::vec3(50.0, 00.0, 50.0);
        }
        initRegion = AABB(glm::vec3(2.0, 0.2, 2.0), glm::vec3(29.7, 0.2 + 32, 29.7));
        for (int i = 0; i < halfN; ++i) {
            int dim[3] = {
                i % dimParticlesCnt[0],
                i / dimParticlesCnt[0] % 32,
                i / (dimParticlesCnt[0] * 32)
            };
            position[i + halfN] = initRegion.pMin + glm::vec3(
                dim[0] * dimGap[0], dim[1] * dimGap[1], dim[2] * dimGap[2]
            ) + 0.2f * glm::vec3(float(rand())/RAND_MAX, float(rand())/RAND_MAX, float(rand())/RAND_MAX);
            velocity[i + halfN] = glm::vec3(-50.0, 0.0, -50.0);
        }
    }
    else if (fluidSceneIdx == 4) {
        initRegion = AABB(glm::vec3(-24.0, 0.2, -16.0), glm::vec3(8.0, 0.2 + 32, 16.0));
        dimGap[0] = (initRegion.pMax.x - initRegion.pMin.x) / dimParticlesCnt[0];
        dimGap[1] = (initRegion.pMax.y - initRegion.pMin.y) / dimParticlesCnt[1];
        dimGap[2] = (initRegion.pMax.z - initRegion.pMin.z) / dimParticlesCnt[2];
        for (int i = 0; i < para.nParticles; ++i) {
            int dim[3] = {
                i % dimParticlesCnt[0],
                i / dimParticlesCnt[0] % dimParticlesCnt[1],
                i / (dimParticlesCnt[0] * dimParticlesCnt[1])
            };
            position[i] = initRegion.pMin + glm::vec3(
                dim[0] * dimGap[0], dim[1] * dimGap[1], dim[2] * dimGap[2]
            ) + 0.2f * glm::vec3(float(rand())/RAND_MAX, float(rand())/RAND_MAX, float(rand())/RAND_MAX);
            velocity[i] = glm::vec3(-10.0, 15.0, 0.0);
        }
    }
    else if (fluidSceneIdx == 5) {
        initRegion = AABB(glm::vec3(-16.0, 0.2, -16.0), glm::vec3(16, 0.2 + 32, 16.0));
        dimGap[0] = (initRegion.pMax.x - initRegion.pMin.x) / dimParticlesCnt[0];
        dimGap[1] = (initRegion.pMax.y - initRegion.pMin.y) / dimParticlesCnt[1];
        dimGap[2] = (initRegion.pMax.z - initRegion.pMin.z) / dimParticlesCnt[2];
        for (int i = 0; i < para.nParticles; ++i) {
            int dim[3] = {
                i % dimParticlesCnt[0],
                i / dimParticlesCnt[0] % dimParticlesCnt[1],
                i / (dimParticlesCnt[0] * dimParticlesCnt[1])
            };
            position[i] = initRegion.pMin + glm::vec3(
                dim[0] * dimGap[0], dim[1] * dimGap[1], dim[2] * dimGap[2]
            ) + 0.2f * glm::vec3(float(rand())/RAND_MAX, float(rand())/RAND_MAX, float(rand())/RAND_MAX);
            velocity[i] = 50.0f * glm::normalize(-glm::vec3(-position[i].z, 0.0, position[i].x));
        }
    }
    SolverData cudaSolverData;
    uint32_t cudaParticleVAO;
    uint32_t cudaParticleVBO;
    glGenVertexArrays(1, &cudaParticleVAO);
    glBindVertexArray(cudaParticleVAO);
    glGenBuffers(1, &cudaParticleVBO);
    glBindBuffer(GL_ARRAY_BUFFER, cudaParticleVBO);
    glBufferData(GL_ARRAY_BUFFER, para.nParticles * 3 * sizeof(float), NULL, GL_DYNAMIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    cudaGraphicsResource* positionBuffer;
    cudaGraphicsGLRegisterBuffer(&positionBuffer, cudaParticleVBO, cudaGraphicsRegisterFlagsNone);
    cudaGraphicsMapResources(1, &positionBuffer);
    cudaGraphicsResourceGetMappedPointer((void**)&cudaSolverData.position, nullptr, positionBuffer);

    initCUDASolverDataBuffer(position, velocity, cudaSolverData, para, box, boundary);

    particleProgram.use();
    particleProgram.setMat4("projection", projection);
    particleProgram.setMat4("view", camera.getViewMatrix());
    particleProgram.setFloat("radius", para.radius);
    particleProgram.setFloat("pointCoeff", WINDOW_HEIGHT / (std::tan(glm::radians(fov/ 2))));
    particleProgram.setVec3("camPos", camera.pos);
    particleProgram.setVec3("pointLight.pos", glm::vec3(0.0, 50.0, 40.0));
    particleProgram.setVec3("pointLight.diffuse", glm::vec3(0.8, 0.8, 0.8));
    particleProgram.setVec3("pointLight.specular", glm::vec3(1.0, 1.0, 1.0));


    float quadPlane[] = {
        -1.0f,  1.0f, 0.0f,  0.0f, 1.0f,
        -1.0f, -1.0f, 0.0f,  0.0f, 0.0f,
        1.0f, -1.0f, 0.0f,  1.0f, 0.0f,

        -1.0f,  1.0f, 0.0f,  0.0f, 1.0f,
        1.0f, -1.0f, 0.0f,  1.0f, 0.0f,
        1.0f,  1.0f, 0.0f,  1.0f, 1.0f
    };
    unsigned int quadPlaneVAO;
    glGenVertexArrays(1, &quadPlaneVAO);
    glBindVertexArray(quadPlaneVAO);
    unsigned int quadPlaneVBO;
    glGenBuffers(1, &quadPlaneVBO);
    glBindBuffer(GL_ARRAY_BUFFER, quadPlaneVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(quadPlane), quadPlane, GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);

    Program depthProg("../shaders/depth.vert", "../shaders/depth.frag");
    Program smoothDepthProg("../shaders/smooth.vert", "../shaders/smooth.frag");
    Program thicknessProg("../shaders/thickness.vert", "../shaders/thickness.frag");
    Program normalProg("../shaders/normal.vert", "../shaders/normal.frag");
    Program fluidShadingProg("../shaders/shadingParticle.vert", "../shaders/shadingParticle.frag");
    Program testProg("../shaders/test.vert", "../shaders/test.frag");
    // depth buffer
    unsigned int depthBuffer;
    glGenRenderbuffers(1, &depthBuffer);
    glBindRenderbuffer(GL_RENDERBUFFER, depthBuffer);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT32F, WINDOW_WIDTH, WINDOW_HEIGHT);
    glBindRenderbuffer(GL_RENDERBUFFER, 0);
    // depthTex
    unsigned int depthTexImg[2];
    unsigned int srcDepTexIdx = 0;
    unsigned int dstDepTexIdx = 1;
    glGenTextures(1, &depthTexImg[0]);
    glBindTexture(GL_TEXTURE_2D, depthTexImg[0]);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, WINDOW_WIDTH, WINDOW_HEIGHT, 0, GL_RED, GL_FLOAT, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glGenTextures(1, &depthTexImg[1]);
    glBindTexture(GL_TEXTURE_2D, depthTexImg[1]);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, WINDOW_WIDTH, WINDOW_HEIGHT, 0, GL_RED, GL_FLOAT, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    // thicknessTex
    unsigned int thicknessTexImg;
    glGenTextures(1, &thicknessTexImg);
    glBindTexture(GL_TEXTURE_2D, thicknessTexImg);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, WINDOW_WIDTH, WINDOW_HEIGHT, 0, GL_RED, GL_FLOAT, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    // normalTex
    unsigned int normalTexImg;
    glGenTextures(1, &normalTexImg);
    glBindTexture(GL_TEXTURE_2D, normalTexImg);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB32F, WINDOW_WIDTH, WINDOW_HEIGHT, 0, GL_RGB, GL_FLOAT, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    // backgroundTex
    unsigned int backgroudTexImg;
    glGenTextures(1, &backgroudTexImg);
    glBindTexture(GL_TEXTURE_2D, backgroudTexImg);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, WINDOW_WIDTH, WINDOW_HEIGHT, 0, GL_RGBA, GL_FLOAT, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    // fluid render FBO
    GLenum depthTexDrawTarget[2] = { GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1 };
    GLenum thicknessTexDrawTarget = GL_COLOR_ATTACHMENT2;
    GLenum normalTexDrawTarget = GL_COLOR_ATTACHMENT3;
    GLenum backgroundTexDrawTaget = GL_COLOR_ATTACHMENT4;
    unsigned int fluidRenderFBO;
    glGenFramebuffers(1, &fluidRenderFBO);
    glBindFramebuffer(GL_FRAMEBUFFER, fluidRenderFBO);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, depthBuffer);
    glFramebufferTexture2D(GL_FRAMEBUFFER, depthTexDrawTarget[0], GL_TEXTURE_2D, depthTexImg[0], 0);
    glFramebufferTexture2D(GL_FRAMEBUFFER, depthTexDrawTarget[1], GL_TEXTURE_2D, depthTexImg[1], 0);
    glFramebufferTexture2D(GL_FRAMEBUFFER, thicknessTexDrawTarget, GL_TEXTURE_2D, thicknessTexImg, 0);
    glFramebufferTexture2D(GL_FRAMEBUFFER, normalTexDrawTarget, GL_TEXTURE_2D, normalTexImg, 0);
    glFramebufferTexture2D(GL_FRAMEBUFFER, backgroundTexDrawTaget, GL_TEXTURE_2D, backgroudTexImg, 0);
    checkFramebufferComplete();
    glBindFramebuffer(GL_FRAMEBUFFER, 0);

    // fluid shader parameter
    float tanHalfFov = std::tan(glm::radians(fov/ 2));
    float pointCoeff = WINDOW_HEIGHT / tanHalfFov;
    float Fy = projection[0][0];
    float Fx = projection[1][1];
    float ZERO = 0.0f;
    float INF = 1000.0f;
    float BLACK[] = {0.0f, 0.0f, 0.0f};
    float smoothingDt = 0.001;

    depthProg.use();
    depthProg.setMat4("projection", projection);
    depthProg.setMat4("view", view);
    depthProg.setFloat("radius", para.radius);
    depthProg.setFloat("pointCoeff", pointCoeff);

    smoothDepthProg.use();
    smoothDepthProg.setFloat("width", (float)WINDOW_WIDTH);
    smoothDepthProg.setFloat("height", (float)WINDOW_HEIGHT);
    smoothDepthProg.setFloat("Fx", Fx);
    smoothDepthProg.setFloat("Fy", Fy);
    smoothDepthProg.setFloat("smoothingDt", smoothingDt);

    thicknessProg.use();
    thicknessProg.setMat4("projection", projection);
    thicknessProg.setMat4("view", view);
    thicknessProg.setFloat("radius", para.radius);
    thicknessProg.setFloat("pointCoeff", pointCoeff);

    normalProg.use();
    normalProg.setFloat("width", (float)WINDOW_WIDTH);
    normalProg.setFloat("height", (float)WINDOW_HEIGHT);
    normalProg.setFloat("Fx", Fx);
    normalProg.setFloat("Fy", Fy);

    fluidShadingProg.use();
    fluidShadingProg.setVec3("pointLight.pos", glm::vec3(0.0, 50.0, -20.0));
    fluidShadingProg.setVec3("pointLight.diffuse", glm::vec3(0.8, 0.8, 0.8));
    fluidShadingProg.setVec3("pointLight.specular", glm::vec3(1.0, 1.0, 1.0));
    fluidShadingProg.setFloat("roi", 1.3f);
    fluidShadingProg.setFloat("width", WINDOW_WIDTH);
    fluidShadingProg.setFloat("height", WINDOW_HEIGHT);
    fluidShadingProg.setFloat("tanHalfFov", tanHalfFov);
    fluidShadingProg.setMat4("viewMat", view);

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_PROGRAM_POINT_SIZE);

    while(!glfwWindowShouldClose(window)) {
        // rendering
        if (fluidRenderMode == 0) {
            auto start = std::chrono::system_clock::now();
            glBindFramebuffer(GL_FRAMEBUFFER, fluidRenderFBO);
            srcDepTexIdx = 0;
            dstDepTexIdx = 1;

            // render depthTexImg
            glEnable(GL_DEPTH_TEST);
            glDisable(GL_BLEND);
            glClear(GL_DEPTH_BUFFER_BIT);
            glClearTexImage(depthTexImg[0], 0, GL_RED, GL_FLOAT, &ZERO);
            glClearTexImage(depthTexImg[1], 0, GL_RED, GL_FLOAT, &ZERO);
            checkGLError();
            glDrawBuffers(1, &depthTexDrawTarget[srcDepTexIdx]);
            depthProg.use();
            glBindVertexArray(cudaParticleVAO);
            glDrawArrays(GL_POINTS, 0, para.nParticles);
            // smooth depthTexImg
            glDisable(GL_DEPTH_TEST);
            for (int iter = 0; iter < 60; ++iter) {
                smoothDepthProg.use();
                glDrawBuffers(1, &depthTexDrawTarget[dstDepTexIdx]);
                glActiveTexture(GL_TEXTURE0);
                glBindTexture(GL_TEXTURE_2D, depthTexImg[srcDepTexIdx]);
                glBindVertexArray(quadPlaneVAO);
                glDrawArrays(GL_TRIANGLES, 0, 6);
                std::swap(srcDepTexIdx, dstDepTexIdx);
            }
            // render thicknessTexImg
            glEnable(GL_BLEND);
            glDisable(GL_DEPTH_TEST);
            glClearTexImage(thicknessTexImg, 0, GL_RED, GL_FLOAT, &ZERO);
            glDrawBuffers(1, &thicknessTexDrawTarget);
            glBlendEquationSeparate(GL_FUNC_ADD, GL_FUNC_ADD);
            glBlendFuncSeparate(GL_ONE, GL_ONE, GL_ONE, GL_ONE);
            thicknessProg.use();
            glBindVertexArray(cudaParticleVAO);
            glDrawArrays(GL_POINTS, 0, para.nParticles);
            // render normalTexImg
            glDisable(GL_DEPTH_TEST);
            glDisable(GL_BLEND);
            glDrawBuffers(1, &normalTexDrawTarget);
            glClearTexImage(normalTexImg, 0, GL_RGB, GL_FLOAT, BLACK);
            normalProg.use();
            glActiveTexture(GL_TEXTURE0);
            glBindTexture(GL_TEXTURE_2D, depthTexImg[srcDepTexIdx]);
            glBindVertexArray(quadPlaneVAO);
            glDrawArrays(GL_TRIANGLES, 0, 6);
            // render backgroudTexImg
            glDrawBuffers(1, &backgroundTexDrawTaget);
            glClearColor(0.75f, 0.75f, 0.75f, 1.0f);
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            glEnable(GL_DEPTH_TEST);
            backgroundProgram.use();
            glBindVertexArray(backgroundVAO);
            glDrawArrays(GL_TRIANGLES, 0, 24);
            // fluid shading
            glBindFramebuffer(GL_FRAMEBUFFER, 0);
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            glEnable(GL_DEPTH_TEST);
            fluidShadingProg.use();
            glActiveTexture(GL_TEXTURE0);
            glBindTexture(GL_TEXTURE_2D, depthTexImg[srcDepTexIdx]);
            glActiveTexture(GL_TEXTURE1);
            glBindTexture(GL_TEXTURE_2D, thicknessTexImg);
            glActiveTexture(GL_TEXTURE2);
            glBindTexture(GL_TEXTURE_2D, normalTexImg);
            glActiveTexture(GL_TEXTURE3);
            glBindTexture(GL_TEXTURE_2D, backgroudTexImg);
            fluidShadingProg.setInt("zValTex", 0);
            fluidShadingProg.setInt("thicknessTex", 1);
            fluidShadingProg.setInt("normalTex", 2);
            fluidShadingProg.setInt("backgroundTex", 3);
            glBindVertexArray(quadPlaneVAO);
            glDrawArrays(GL_TRIANGLES, 0, 6);
            // test
            //glBindFramebuffer(GL_FRAMEBUFFER, 0);
            //glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            ////glEnable(GL_DEPTH_TEST);
            //testProg.use();
            //glActiveTexture(GL_TEXTURE0);
            //glBindTexture(GL_TEXTURE_2D, normalTexImg);
            //glBindVertexArray(quadPlaneVAO);
            //glDrawArrays(GL_TRIANGLES, 0, 6);
        }
        else if (fluidRenderMode== 1) {
            glClearColor(0.5f, 0.5f, 0.5f, 1.0f);
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            glDisable(GL_BLEND);
            glEnable(GL_DEPTH_TEST);
            backgroundProgram.use();
            glBindVertexArray(backgroundVAO);
            glDrawArrays(GL_TRIANGLES, 0, 24);
            particleProgram.use();
            glBindVertexArray(cudaParticleVAO);
            glDrawArrays(GL_POINTS, 0, para.nParticles);

        }
        // simulation
        auto start = std::chrono::system_clock::now();
        cudaPBFSolve(cudaSolverData, para);
        auto end = std::chrono::system_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        std::cout << double(duration.count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den << std::endl;
        glfwSwapBuffers(window);
        glfwPollEvents();

    }//
    cudaGraphicsUnmapResources(1, &positionBuffer);
    //freeCUDASolverDataBuffer(cudaSolverData);
    glDeleteBuffers(1, &backgroundVBO);
    glDeleteVertexArrays(1, &backgroundVAO);
}

void RenderApp::cleanUp() {
    glfwTerminate();
}

void RenderApp::checkFramebufferComplete()
{
	GLenum errorCode = glCheckFramebufferStatus(GL_FRAMEBUFFER);
	const char *errorInfo = NULL;
	switch (errorCode) {
	case GL_FRAMEBUFFER_UNDEFINED: errorInfo = "GL_FRAMEBUFFER_UNDEFINED"; break;
	case GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT: errorInfo = "GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT"; break;
	case GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT: errorInfo = "GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT"; break;
	case GL_FRAMEBUFFER_INCOMPLETE_DRAW_BUFFER: errorInfo = "GL_FRAMEBUFFER_INCOMPLETE_DRAW_BUFFER"; break;
	case GL_FRAMEBUFFER_INCOMPLETE_READ_BUFFER: errorInfo = "GL_FRAMEBUFFER_INCOMPLETE_READ_BUFFER"; break;
	case GL_FRAMEBUFFER_UNSUPPORTED: errorInfo = "GL_FRAMEBUFFER_UNSUPPORTED"; break;
	case GL_FRAMEBUFFER_INCOMPLETE_MULTISAMPLE: errorInfo = "GL_FRAMEBUFFER_INCOMPLETE_MULTISAMPLE"; break;
	case GL_FRAMEBUFFER_INCOMPLETE_LAYER_TARGETS: errorInfo = "GL_FRAMEBUFFER_INCOMPLETE_LAYER_TARGETS"; break;
	}
	if (errorInfo) {
		fprintf(stderr, "OpenGL Framebuffer Error #%d: %s\n", errorCode, errorInfo);
		exit(-1);
	}
	else {
		printf("Framebuffer complete check ok\n");
	}
}

void RenderApp::checkGLError() {
    GLenum errorCode;
	const char *errorInfo;
	if ((errorCode = glGetError()) != GL_NO_ERROR) {
		switch (errorCode) {
			case GL_INVALID_OPERATION:      errorInfo = "INVALID_OPERATION";      break;
			case GL_INVALID_ENUM:           errorInfo = "INVALID_ENUM";           break;
			case GL_INVALID_VALUE:          errorInfo = "INVALID_VALUE";          break;
			case GL_OUT_OF_MEMORY:          errorInfo = "OUT_OF_MEMORY";          break;
			case GL_INVALID_FRAMEBUFFER_OPERATION:  errorInfo = "INVALID_FRAMEBUFFER_OPERATION";  break;
		}
		fprintf(stderr, "OpenGL Error #%d: %s\n", errorCode, errorInfo);
		exit(-1);
	}
}