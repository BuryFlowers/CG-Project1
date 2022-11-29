#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include "camera.h"
#include "shader.h"

using namespace glm;

//frame time
float deltaTime = 0.0f;
float lastFrame = 0.0f;

//camera
camera* cam;
float camSpeed = 5.0f;

//cursor
vec2 lastCursor;
float cursorSensitivity = 0.1f;
bool firstMouse = true;
float blend = 0.2f;

void framebuffer_size_callback(GLFWwindow* window, int width, int height) {

	glViewport(0, 0, width, height);

}

void mouse_callback(GLFWwindow* window, double posX, double posY) {

	float xpos = static_cast<float>(posX);
	float ypos = static_cast<float>(posY);

	if (firstMouse)
	{
		lastCursor = vec2(xpos, ypos);
		firstMouse = false;
	}

	float xoffset = xpos - lastCursor.x;
	float yoffset = lastCursor.y - ypos; // reversed since y-coordinates go from bottom to top

	lastCursor.x = xpos;
	lastCursor.y = ypos;

	cam->ProcessMouseMovement(xoffset, yoffset, true);

}

void scroll_callback(GLFWwindow* window, double offsetX, double offsetY) {

	cam->ProcessMouseScroll(offsetY);

}

void processInput(GLFWwindow* window) {

	deltaTime = glfwGetTime() - lastFrame;
	lastFrame = glfwGetTime();

	if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
		glfwSetWindowShouldClose(window, true);

	if (glfwGetKey(window, GLFW_KEY_UP) == GLFW_PRESS) {

		if (blend < 0.999f) blend += 0.001f;

	}

	if (glfwGetKey(window, GLFW_KEY_DOWN) == GLFW_PRESS) {

		if (blend > 0.001f) blend -= 0.001f;

	}

	if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS) cam->ProcessKeyboard(FORWARD, deltaTime);
	if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS) cam->ProcessKeyboard(BACKWARD, deltaTime);
	if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS) cam->ProcessKeyboard(LEFT, deltaTime);
	if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS) cam->ProcessKeyboard(RIGHT, deltaTime);


}