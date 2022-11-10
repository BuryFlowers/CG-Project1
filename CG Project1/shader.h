#pragma once
#ifndef SHADER
#define SHADER

#include <glad/glad.h>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

using namespace std;

class Shader {

public:

	unsigned int ID;

	//constructor read shaders and construct them
	Shader(const char* vertexShaderPath, const char* fragmentShaderPath);
	//use/activate program
	void use();
	//uniform tools
	void setBool(const string& name, bool value) const;
	void setInt(const string& name, int value) const;
	void setFloat(const string& name, float value) const;
	void setMat4(const string& name, glm::mat4 value) const;
	void setVec3(const string& name, float x, float y, float z) const;
	void setVec3(const string& name, glm::vec3 value) const;

private:


};

#endif // !SHADER
