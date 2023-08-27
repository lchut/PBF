cmake_minimum_required(VERSION 3.10)

project(PBF LANGUAGES CXX C CUDA)

set(SOURCES
    ${PROJECT_SOURCE_DIR}/src/main.cpp
    ${PROJECT_SOURCE_DIR}/src/renderapp.cpp
    ${PROJECT_SOURCE_DIR}/src/program.cpp
    ${PROJECT_SOURCE_DIR}/src/AABB.cpp
    ${PROJECT_SOURCE_DIR}/src/solver.cpp
    ${PROJECT_SOURCE_DIR}/src/glad.c
)

set(SOURCES_CUDA
    ${PROJECT_SOURCE_DIR}/kernels/cudaPBF.cu
    ${PROJECT_SOURCE_DIR}/kernels/cudaPBFSolver.cu
)



add_executable(PBF  ${SOURCES} ${SOURCES_CUDA})

find_package(CUDA REQUIRED)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(PBF OpenMP::OpenMP_CXX)
else()
    message("Could not find OpenMp")
endif()

target_compile_options(PBF PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:
    -Xcudafe
    --diag_suppress=esa_on_defaulted_function_ignored
    >
)

target_include_directories(PBF PUBLIC
    ${CUDA_INCLUDE_DIR}/
    ${PROJECT_SOURCE_DIR}/include/
    ${PROJECT_SOURCE_DIR}/kernels/
    ${PROJECT_SOURCE_DIR}/ext/
)

set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
add_subdirectory(./ext/glfw)
target_link_libraries(PBF glfw)

find_package(OpenGL REQUIRED)
target_link_libraries(PBF OpenGL::GL)