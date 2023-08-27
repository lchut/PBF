#include "renderapp.h"

int main(int argc, char* argv[]) {
    RenderApp app;
    try {
        int sceneIdx = 4;
        int renderMode = 0;
        if (argc == 2) {
            sceneIdx = std::stoi(argv[1]);
        }
        else if (argc == 3) {
            sceneIdx = std::stoi(argv[1]);
            renderMode = std::stoi(argv[2]);
        }
        if (sceneIdx < 1 || sceneIdx > 5) {
            sceneIdx = 1;
        }
        if (renderMode < 0 || renderMode > 1) {
            renderMode = 0;
        }
        fluidSceneIdx = sceneIdx;
        fluidRenderMode = renderMode;
        app.run();
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}