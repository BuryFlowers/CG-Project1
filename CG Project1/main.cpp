#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include "camera.h"
#include "shader.h"
#include "tiny_gltf.h"

#define WIDTH 1280
#define HEIGHT 720

using namespace glm;
using namespace tinygltf;
using namespace std;

//frame time
extern float deltaTime;
extern float lastFrame;

//camera
extern camera* cam;
extern float camSpeed;

//cursor
extern vec2 lastCursor;
extern float cursorSensitivity;
extern bool firstMouse;
extern float blend;

void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void mouse_callback(GLFWwindow* window, double posX, double posY);
void scroll_callback(GLFWwindow* window, double offsetX, double offsetY);
void processInput(GLFWwindow* window);

int main() {

	Model model;
	TinyGLTF loader;
	string err;
	string warn;

	bool ret = loader.LoadASCIIFromFile(&model, &err, &warn, "models/blue_archivekasumizawa_miyu/scene.gltf");

	if (!warn.empty()) {
		printf("Warn: %s\n", warn.c_str());
	}

	if (!err.empty()) {
		printf("Err: %s\n", err.c_str());
	}

	if (!ret) {
		printf("Failed to parse glTF\n");
		return -1;
	}

	printf("Model is loaded!\n");

	//camera set up
	cam = new camera(vec3(0.0f, 0.0f, 10.0f), vec3(0.0f, 1.0f, 0.0f), -90.0f, 0.0f);
	cam->MouseSensitivity = cursorSensitivity;
	cam->MovementSpeed = camSpeed;
	lastCursor = vec2(WIDTH * 0.5f, HEIGHT * 0.5f);

	//get view transform matrix
	mat4 view = cam->GetViewMatrix();

	glfwInit();
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	GLFWwindow* window = glfwCreateWindow(WIDTH, HEIGHT, "CG Project1", NULL, NULL);

	if (window == NULL) {

		printf("Failed to create a window\n");
		glfwTerminate();
		return -1;

	}
	glfwMakeContextCurrent(window);

	if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {

		printf("Failed to initialize GLAD\n");
		return -1;

	}

	glViewport(0, 0, WIDTH, HEIGHT);
	glfwSetWindowSizeCallback(window, framebuffer_size_callback);
	glfwSetCursorPosCallback(window, mouse_callback);
	glfwSetScrollCallback(window, scroll_callback);
	glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

	Shader shader("", "");

	while (!glfwWindowShouldClose(window)) {

		glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		processInput(window);

		glfwSwapBuffers(window);
		glfwPollEvents();

	}

	glfwTerminate();
	return 0;

}