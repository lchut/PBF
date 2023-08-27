#ifndef PBF_RENDERAPP_H_
#define PBF_RENDERAPP_H_

#include "config.h"

#include <iostream>
#include <stdexcept>
#include <string>

const uint32_t WINDOW_WIDTH = 800;
const uint32_t WINDOW_HEIGHT = 800;

extern int fluidSceneIdx;
extern int fluidRenderMode;

class RenderApp {
public:
    RenderApp() : window(nullptr) {}
    void run();
private:
    void initWindow();
    void initOpenGL();
    void mainLoop();
    void cleanUp();
    void checkFramebufferComplete();
    void checkGLError();
    GLFWwindow* window;
};
#endif