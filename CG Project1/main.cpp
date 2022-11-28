#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include "cuda_runtime_api.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "camera.h"
#include "shader.h"
#include "triangle.h"

#include <omp.h>
#define NOMINMAX
#include <windows.h>  
# define	TIMING_BEGIN \
	{double tmp_timing_start = omp_get_wtime();

# define	TIMING_END(message) \
	{double tmp_timing_finish = omp_get_wtime();\
	double  tmp_timing_duration = tmp_timing_finish - tmp_timing_start;\
	printf("\r%s: %2.5f ms           ", (message), tmp_timing_duration * 1000);}}


#define WIDTH 1280
#define HEIGHT 720

using namespace glm;
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

//mesh
extern vector<vec3> position_queue;
extern vector<vec3> normal_queue;
extern vector<vec2> uv_queue;
extern vector<triangle> indices_queue;

triangles mesh;

extern void load_obj_file(const char* filename);

//cuda
extern void cuda_init(float* h_transform, triangles mesh);

//opengl callback
void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void mouse_callback(GLFWwindow* window, double posX, double posY);
void scroll_callback(GLFWwindow* window, double offsetX, double offsetY);
void processInput(GLFWwindow* window);

int main() {

	load_obj_file("data/Ornat_czarny_03_FIN.obj");
	mesh.position = new float[position_queue.size() * 3];
	mesh.normal = new float[normal_queue.size() * 3];
	mesh.uv = new float[uv_queue.size() * 2];
	mesh.triangle_indices = new triangle[uv_queue.size()];

	for (int i = 0; i < position_queue.size(); i++) {

		mesh.position[i * 3 + 0] = position_queue[i].x;
		mesh.position[i * 3 + 1] = position_queue[i].y;
		mesh.position[i * 3 + 2] = position_queue[i].z;
		mesh.normal[i * 3 + 0] = normal_queue[i].x;
		mesh.normal[i * 3 + 1] = normal_queue[i].y;
		mesh.normal[i * 3 + 2] = normal_queue[i].z;
		mesh.uv[i * 2 + 0] = uv_queue[i].x;
		mesh.uv[i * 2 + 1] = uv_queue[i].y;
		mesh.triangle_indices[i] = indices_queue[i];

	}

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
	glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);

	//Shader shader("", "");

	while (!glfwWindowShouldClose(window)) {

		TIMING_BEGIN
		glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		processInput(window);

		mat4 transform = perspective(radians(cam->Zoom), (float)WIDTH / (float)HEIGHT, 0.1f, 100.0f) * cam->GetViewMatrix();
		float* h_transform = new float[16];
		for (int i = 0; i < 16; i++) h_transform[i] = transform[i / 4][i % 4];
		
		cuda_init(h_transform, mesh);


		glfwSwapBuffers(window);
		glfwPollEvents();

		TIMING_END("Time per frame")
	}

	glfwTerminate();
	return 0;

}