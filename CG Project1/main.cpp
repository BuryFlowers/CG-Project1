#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <iostream>
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

unsigned char* PixelBuffer = new unsigned char[WIDTH * HEIGHT * 3];

float vertices[] = {
	//  ---- position ----    ---- color ----     - texcoord -
		 1.0f,  1.0f, 0.0f,   1.0f, 0.0f, 0.0f,   0.0f, 0.0f,   
		 1.0f, -1.0f, 0.0f,   0.0f, 1.0f, 0.0f,   0.0f, 1.0f,   
		-1.0f, -1.0f, 0.0f,   0.0f, 0.0f, 1.0f,   1.0f, 1.0f,   
		-1.0f,  1.0f, 0.0f,   1.0f, 1.0f, 0.0f,   1.0f, 0.0f   
};
unsigned int indices[] = {
		0, 1, 3, // first triangle
		1, 2, 3  // second triangle
};

extern void load_obj_file(const char* filename);

//cuda
extern void cuda_init(float* h_transform, triangles mesh);

//opengl callback
void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void mouse_callback(GLFWwindow* window, double posX, double posY);
void scroll_callback(GLFWwindow* window, double offsetX, double offsetY);
void processInput(GLFWwindow* window);

int main() {

	load_obj_file("data/cube.obj");

	mesh.position = new float[position_queue.size() * 3];
	for (int i = 0; i < position_queue.size(); i++) {

		mesh.position[i * 3] = position_queue[i].x;
		mesh.position[i * 3 + 1] = position_queue[i].y;
		mesh.position[i * 3 + 2] = position_queue[i].z;
		//printf("%.2lf %.2lf %.2lf\n", mesh.position[i * 3], mesh.position[i * 3 + 1], mesh.position[i * 3 + 2]);

	}

	mesh.normal = new float[normal_queue.size() * 3];
	for (int i = 0; i < normal_queue.size(); i++) {

		mesh.normal[i * 3] = normal_queue[i].x;
		mesh.normal[i * 3 + 1] = normal_queue[i].y;
		mesh.normal[i * 3 + 2] = normal_queue[i].z;

	}

	mesh.uv = new float[uv_queue.size() * 2];
	for (int i = 0; i < uv_queue.size(); i++) {

		mesh.uv[i * 2] = uv_queue[i].x;
		mesh.uv[i * 2 + 1] = uv_queue[i].y;

	}

	mesh.triangle_indices = new triangle[indices_queue.size()];
	for (int i = 0; i < indices_queue.size(); i++)
		mesh.triangle_indices[i] = indices_queue[i];

	//camera set up
	cam = new camera(vec3(0.0f, 0.0f, -10.0f), vec3(0.0f, 1.0f, 0.0f), -90.0f, 0.0f);
	cam->MouseSensitivity = cursorSensitivity;
	cam->MovementSpeed = camSpeed;
	lastCursor = vec2(WIDTH * 0.5f, HEIGHT * 0.5f);

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

	Shader shader("default.vs", "default.fs");

	// init VBO VAO EBO
	unsigned int VBO, VAO, EBO;
	glGenVertexArrays(1, &VAO);
	glGenBuffers(1, &VBO);
	glGenBuffers(1, &EBO);

	glBindVertexArray(VAO);

	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);

	// position attribute
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(0);
	// color attribute
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(3 * sizeof(float)));
	glEnableVertexAttribArray(1);
	// texture coord attribute
	glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(6 * sizeof(float)));
	glEnableVertexAttribArray(2);

	//create texture canvas
	unsigned int texture;
	glGenTextures(1, &texture);
	glBindTexture(GL_TEXTURE_2D, texture);
	// set the texture wrapping parameters
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);	// set texture wrapping to GL_REPEAT (default wrapping method)
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	// set texture filtering parameters
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

	// set sampler to sample texture0 
	shader.use();
	glActiveTexture(GL_TEXTURE0);
	shader.setInt("texture1", 0);

	glDisable(GL_DEPTH_TEST);

	while (!glfwWindowShouldClose(window)) {

		TIMING_BEGIN

		processInput(window);
		
		glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT);

		for (int i = 0; i < WIDTH * HEIGHT * 3; i++) PixelBuffer[i] = 0;
		mat4 transform = perspective(radians(cam->Zoom), (float)WIDTH / (float)HEIGHT, 0.1f, 100.0f) * cam->GetViewMatrix();
		for (int i = 0; i < indices_queue.size(); i++) {

			int A = mesh.triangle_indices[i].A();
			int B = mesh.triangle_indices[i].B();
			int C = mesh.triangle_indices[i].C();
			float* position = mesh.position;
			vec4 p[3];
			p[0] = vec4(position[A * 3], position[A * 3 + 1], position[A * 3 + 2], 1.0f);
			p[1] = vec4(position[B * 3], position[B * 3 + 1], position[B * 3 + 2], 1.0f);
			p[2] = vec4(position[C * 3], position[C * 3 + 1], position[C * 3 + 2], 1.0f);
			p[0] = transform * p[0];
			p[1] = transform * p[1];
			p[2] = transform * p[2];

			for (int j = 0; j < 3; j++) {

				vec3 t = vec3((p[(j + 1) % 3].x + 1) / 2.0 * WIDTH, (p[(j + 1) % 3].y + 1) / 2.0 * HEIGHT, p[(j + 1) % 3].z);
				vec3 s = vec3((p[j].x + 1) / 2.0 * WIDTH, (p[j].y + 1) / 2.0 * HEIGHT, p[j].z);
				if (t.x < s.x) swap(s, t);
				float k = (t.y - s.y) / (t.x - s.x);

				if (fabs(k) > 1.0f) {

					for (int y = std::max((int)s.y, 0); y < std::min((int)t.y, HEIGHT); y++) {

						int x = ((float)y - s.y) / k + s.x;
						float z = ((float)y - s.y) / (t.y - s.y) * (t.z - s.z) + s.z;
						int index = (int)(x + y * WIDTH);
						if (index < WIDTH * HEIGHT && index > 0 && z > -1) {

							PixelBuffer[index * 3] = 255;
							PixelBuffer[index * 3 + 1] = 255;
							PixelBuffer[index * 3 + 2] = 255;

						}

					}

				}

				else {

					for (int x = std::max((int)s.x, 0); x < std::min((int)t.x, WIDTH); x++) {

						int y = ((float)x - s.x) * k + s.y;
						float z = ((float)x - s.x) / (t.x - s.x) * (t.z - s.z) + s.z;
						int index = (int)(x + y * WIDTH);
						if (index < WIDTH * HEIGHT && index > 0 && z > -1) {

							PixelBuffer[index * 3] = 255;
							PixelBuffer[index * 3 + 1] = 255;
							PixelBuffer[index * 3 + 2] = 255;

						}

					}

				}

			}

		}
		
		//cuda_init(h_transform, mesh);

		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, WIDTH, HEIGHT, 0, GL_RGB, GL_UNSIGNED_BYTE, PixelBuffer);

		// render canvas
		shader.use();
		glBindVertexArray(VAO);
		glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

		glfwSwapBuffers(window);
		glfwPollEvents();

		TIMING_END("Time per frame")
	}

	glfwTerminate();
	return 0;

}