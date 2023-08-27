#include "program.h"

void checkCompileErrors(GLuint shader, std::string type)
{
    GLint success;
    GLchar infoLog[1024];
    if(type != "PROGRAM")
    {
        glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
        if(!success)
        {
            glGetShaderInfoLog(shader, 1024, NULL, infoLog);
            std::cout << "ERROR::SHADER_COMPILATION_ERROR of type: " << type << "\n" << infoLog << "\n -- --------------------------------------------------- -- " << std::endl;
        }
    }
    else
    {
        glGetProgramiv(shader, GL_LINK_STATUS, &success);
        if(!success)
        {
            glGetProgramInfoLog(shader, 1024, NULL, infoLog);
            std::cout << "ERROR::PROGRAM_LINKING_ERROR of type: " << type << "\n" << infoLog << "\n -- --------------------------------------------------- -- " << std::endl;
        }
    }
}

Program::Program(const std::string& vertexShaderFileName, const std::string& fragmentShaderFileName) {
    uint32_t vertexShaderId = createShader(vertexShaderFileName, GL_VERTEX_SHADER);
    uint32_t fragmentShaderId = createShader(fragmentShaderFileName, GL_FRAGMENT_SHADER);
    id = glCreateProgram();
    glAttachShader(id, vertexShaderId);
    glAttachShader(id, fragmentShaderId);
    glLinkProgram(id);
    checkCompileErrors(id, "PROGRAM");
    glDeleteShader(vertexShaderId);
    glDeleteShader(fragmentShaderId);
}

std::string Program::readFile(const std::string& filename) const {
    std::ifstream file;
    file.exceptions(std::ifstream::failbit | std::ifstream::badbit);
    std::string ret;
    try {
        file.open(filename);
        std::stringstream fileStream;
        fileStream << file.rdbuf();
        file.close();
        ret = fileStream.str();
    }
    catch(std::ifstream::failure& e) {
        std::cout << "ERROR::SHADER::FILE_NOT_SUCCESFULLY_READ: " << e.what() << std::endl;
    }
    return ret;
}

uint32_t Program::createShader(const std::string& shaderFileName, GLenum shaderType) const {
    std::string code = readFile(shaderFileName);
    const char* shaderCode = code.c_str();
    uint32_t shaderId = glCreateShader(shaderType);
    glShaderSource(shaderId, 1, &shaderCode, nullptr);
    glCompileShader(shaderId);
    checkCompileErrors(shaderId, "VERTEX");
    return shaderId;
}