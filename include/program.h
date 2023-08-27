#ifndef PBF_PROGRAM_H_
#define PBF_PROGRAM_H_

#include "config.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <cstring>
#include <string>

class Program {
public:
    Program(const std::string& vertexShaderFileName, const std::string& fragmentShaderFileName);
    ~Program() {
        glDeleteProgram(id);
    }
    
    inline void use() { glUseProgram(id); }

    inline void setBool(const std::string &name, bool value) const
    {         
        glUniform1i(glGetUniformLocation(id, name.c_str()), (int)value); 
    }
    
    inline void setInt(const std::string &name, int value) const
    { 
        glUniform1i(glGetUniformLocation(id, name.c_str()), value); 
    }

    inline void setFloat(const std::string &name, float value) const
    { 
        glUniform1f(glGetUniformLocation(id, name.c_str()), value); 
    }

    inline void setVec2(const std::string &name, const glm::vec2 &value) const
    { 
        glUniform2fv(glGetUniformLocation(id, name.c_str()), 1, &value[0]); 
    }
    inline void setVec2(const std::string &name, float x, float y) const
    { 
        glUniform2f(glGetUniformLocation(id, name.c_str()), x, y); 
    }
    
    inline void setVec3(const std::string &name, const glm::vec3 &value) const
    { 
        glUniform3fv(glGetUniformLocation(id, name.c_str()), 1, &value[0]); 
    }

    inline void setVec3(const std::string &name, float x, float y, float z) const
    { 
        glUniform3f(glGetUniformLocation(id, name.c_str()), x, y, z); 
    }

    inline void setVec4(const std::string &name, const glm::vec4 &value) const
    { 
        glUniform4fv(glGetUniformLocation(id, name.c_str()), 1, &value[0]); 
    }
    inline void setVec4(const std::string &name, float x, float y, float z, float w) 
    { 
        glUniform4f(glGetUniformLocation(id, name.c_str()), x, y, z, w); 
    }

    inline void setMat2(const std::string &name, const glm::mat2 &mat) const
    {
        glUniformMatrix2fv(glGetUniformLocation(id, name.c_str()), 1, GL_FALSE, &mat[0][0]);
    }

    inline void setMat3(const std::string &name, const glm::mat3 &mat) const
    {
        glUniformMatrix3fv(glGetUniformLocation(id, name.c_str()), 1, GL_FALSE, &mat[0][0]);
    }

    inline void setMat4(const std::string &name, const glm::mat4 &mat) const
    {
        glUniformMatrix4fv(glGetUniformLocation(id, name.c_str()), 1, GL_FALSE, &mat[0][0]);
    }
private:
    std::string readFile(const std::string& filename) const;
    uint32_t createShader(const std::string& shaderFileName, GLenum shaderType) const;
    uint32_t id;
};

#endif