#include "shader.h"

using namespace std;

Shader::Shader(const char* vertexShaderPath, const char* fragmentShaderPath) {

	string vertexCode;
	string fragmentCode;
	ifstream vShaderFile;
	ifstream fShaderFile;

	//error handle
	vShaderFile.exceptions(ifstream::failbit | ifstream::badbit);
	fShaderFile.exceptions(ifstream::failbit | ifstream::badbit);

	try {

		//open 2 files
		vShaderFile.open(vertexShaderPath);
		fShaderFile.open(fragmentShaderPath);

		stringstream vShaderStream, fShaderStream;
		//read buff to streams
		vShaderStream << vShaderFile.rdbuf();
		fShaderStream << fShaderFile.rdbuf();

		//close 2 filef
		vShaderFile.close();
		fShaderFile.close();

		//turn streams to strings
		vertexCode = vShaderStream.str();
		fragmentCode = fShaderStream.str();

	}

	catch (ifstream::failure error) {

		printf("ERROR::SHADER::FILE_NOT_SUCCESFULLY_READ\n");

	}

	//turn string to char
	const char* vertexShaderCode = vertexCode.c_str();
	const char* fragmentShaderCode = fragmentCode.c_str();

	unsigned int vertex, fragment;
	int success;
	char infoLog[512];

	//create vertex shader
	vertex = glCreateShader(GL_VERTEX_SHADER);
	glShaderSource(vertex, 1, &vertexShaderCode, NULL);
	glCompileShader(vertex);

	//check compile result
	glGetShaderiv(vertex, GL_COMPILE_STATUS, &success);
	if (!success) {

		glGetShaderInfoLog(vertex, 512, NULL, infoLog);
		printf("ERROR::SHADER::VERTEX::COMPILATION_FAILED\n%s\n", infoLog);

	}

	//create fragment shader
	fragment = glCreateShader(GL_FRAGMENT_SHADER);
	glShaderSource(fragment, 1, &fragmentShaderCode, NULL);
	glCompileShader(fragment);

	//check compile result
	glGetShaderiv(fragment, GL_COMPILE_STATUS, &success);
	if (!success) {

		glGetShaderInfoLog(fragment, 512, NULL, infoLog);
		printf("ERROR::SHADER::FRAGMENT::COMPILATION_FAILED\n%s\n", infoLog);

	}

	ID = glCreateProgram();
	glAttachShader(ID, vertex);
	glAttachShader(ID, fragment);
	glLinkProgram(ID);

	//check link result
	glGetShaderiv(ID, GL_LINK_STATUS, &success);
	if (!success) {

		glGetShaderInfoLog(ID, 512, NULL, infoLog);
		printf("ERROR::SHADER::PROGRAM::LINK_FAILED\n%s\n", infoLog);

	}

	glDeleteShader(vertex);
	glDeleteShader(fragment);

}

void Shader::use() {

	glUseProgram(ID);

}

void Shader::setBool(const string& name, bool value) const {

	glUniform1i(glGetUniformLocation(ID, name.c_str()), (int)value);

}

void Shader::setInt(const string& name, int value) const {

	glUniform1i(glGetUniformLocation(ID, name.c_str()), value);

}

void Shader::setFloat(const string& name, float value) const {

	glUniform1f(glGetUniformLocation(ID, name.c_str()), value);

}

void Shader::setMat4(const string& name, glm::mat4 value) const {

	glUniformMatrix4fv(glGetUniformLocation(ID, name.c_str()), 1, GL_FALSE, glm::value_ptr(value));

}

void Shader::setVec3(const string& name, float x, float y, float z) const {

	glUniform3f(glGetUniformLocation(ID, name.c_str()), x, y, z);

}

void Shader::setVec3(const string& name, glm::vec3 value) const {

	glUniform3f(glGetUniformLocation(ID, name.c_str()), value.x, value.y, value.z);

}