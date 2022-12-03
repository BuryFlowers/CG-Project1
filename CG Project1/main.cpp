#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include "cuda_runtime_api.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "camera.h"
#include "shader.h"
#include "triangle.h"
#include "scanline_zbuffer.h"

#include <omp.h>
#define NOMINMAX
#include <windows.h>  
# define	TIMING_BEGIN \
	{double tmp_timing_start = omp_get_wtime();

# define	TIMING_END(message) \
	{double tmp_timing_finish = omp_get_wtime();\
	double  tmp_timing_duration = tmp_timing_finish - tmp_timing_start;\
	printf("\r%s: %2.5f ms           ", (message), tmp_timing_duration * 1000);}}


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
polygon_table* PT[HEIGHT];
edge_table* ET[HEIGHT];
active_polygon_table* APT;
active_edge_table* AET;
float* zbuffer;

unsigned char* PixelBuffer = new unsigned char[WIDTH * HEIGHT * 3];

float vertices[] = {
	//  ---- position ----    ---- color ----     - texcoord -
		 1.0f,  1.0f, 0.0f,   1.0f, 0.0f, 0.0f,   1.0f, 1.0f,   
		 1.0f, -1.0f, 0.0f,   0.0f, 1.0f, 0.0f,   1.0f, 0.0f,   
		-1.0f, -1.0f, 0.0f,   0.0f, 0.0f, 1.0f,   0.0f, 0.0f,   
		-1.0f,  1.0f, 0.0f,   1.0f, 1.0f, 0.0f,   0.0f, 1.0f   
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

inline int getYCoord(float y) { return (int)((y + 1.0f) * 0.5f * HEIGHT + 0.5f); }
inline float getXCoord(float x) { return (x + 1.0f) * 0.5f * WIDTH; }
void attachToEdge(vec3& p1, vec3& p2);
bool shouldBeThrown(vec3 s, vec3 t);

int main() {

	load_obj_file("data/diamond.obj");

	mesh.position = new vec3[position_queue.size()];
	mesh.ndc_position = new vec3[position_queue.size()];
	for (int i = 0; i < position_queue.size(); i++) mesh.position[i] = position_queue[i];

	mesh.normal = new vec3[normal_queue.size()];
	for (int i = 0; i < normal_queue.size(); i++) mesh.normal[i] = normal_queue[i];


	mesh.uv = new vec2[uv_queue.size()];
	for (int i = 0; i < uv_queue.size(); i++) mesh.uv[i] = uv_queue[i];

	mesh.triangle_indices = new triangle[indices_queue.size()];
	for (int i = 0; i < indices_queue.size(); i++) mesh.triangle_indices[i] = indices_queue[i];

	std::srand(std::time(0));
	vec3* color = new vec3[indices_queue.size()];
	for (int i = 0; i < indices_queue.size(); i++) {

		color[i].x = std::rand() % 255;
		color[i].y = std::rand() % 255;
		color[i].z = std::rand() % 255;

	}

	zbuffer = new float[WIDTH];

	//camera set up
	//cam = new camera(vec3(0.0f, 0.0f, -5.0f), vec3(0.0f, 1.0f, 0.0f), 90.0f, 0.0f);
	cam = new camera(vec3(0.738162f, 0.518383f, -1.616352f), vec3(0.0f, 1.0f, 0.0f), 746.304016, 39.900074f);
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
		
		glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT);

		for (int i = 0; i < WIDTH * HEIGHT; i++) {

			PixelBuffer[i * 3] = (int)(255);
			PixelBuffer[i * 3 + 1] = (int)(255);
			PixelBuffer[i * 3 + 2] = (int)(255);

		}
		for (int i = 0; i < HEIGHT; i++) {

			PT[i] = NULL;
			ET[i] = NULL;

		}
		mat4 transform = perspective(radians(cam->Zoom), RATIO, 0.1f, 100.0f) * cam->GetViewMatrix();
		for (int i = 0; i < position_queue.size(); i++) {

			vec4 tmp;
			tmp = transform * vec4(mesh.position[i], 1.0f);
			mesh.ndc_position[i] = vec3(tmp / fabs(tmp.w));

		}

		polygon_table* polygon = new polygon_table[indices_queue.size()];
		edge_table* edge[3];
		edge[0] = new edge_table[indices_queue.size()];
		edge[1] = new edge_table[indices_queue.size()];
		edge[2] = new edge_table[indices_queue.size()];
		for (int i = 0; i < indices_queue.size(); i++) {

			int e[3];
			e[0] = mesh.triangle_indices[i].A();
			e[1] = mesh.triangle_indices[i].B();
			e[2] = mesh.triangle_indices[i].C();
			vec3* ndc_p = mesh.ndc_position;
			vec3 AB = vec3(ndc_p[e[0]] - ndc_p[e[1]]);
			vec3 CB = vec3(ndc_p[e[2]] - ndc_p[e[1]]);
			vec3 normal = cross(AB, CB);
			normal = normalize(normal);
			polygon[i].plane.x = normal.x;
			polygon[i].plane.y = normal.y;
			polygon[i].plane.z = normal.z;
			polygon[i].plane.w = -dot(normal, ndc_p[e[1]]);
			polygon[i].color.x = color[i].x;
			polygon[i].color.y = color[i].y;
			polygon[i].color.z = color[i].z;
			if (fabs(polygon[i].plane.z) < EPSILON) continue;

			int ymax = -1, ymin = HEIGHT;
			for (int j = 0; j < 3; j++) {

				edge[j][i].polygon_id = i;
				vec3 s = ndc_p[e[j]], t = ndc_p[e[(j + 1) % 3]];
				if (s.y > t.y) swap(s, t);
				edge[j][i].out = false;
				if (shouldBeThrown(s, t)) {

					edge[j][i].out = true;
					continue;

				}
				edge[j][i].dx = -RATIO * (s.x - t.x) / (s.y - t.y);
				if (s.x < t.x) attachToEdge(s, t);
				else attachToEdge(t, s);
				int top = getYCoord(t.y), down = getYCoord(s.y);

				edge[j][i].dy = top - down;
				edge[j][i].x = getXCoord(t.x);
				edge[j][i].z = t.z;
				if (edge[j][i].dy == 0) {

					edge[j][i].out = true;
					continue;

				}
				if (top > ymax) ymax = top;
 				if (down < ymin) ymin = down;
				edge[j][i].next = NULL;

				if (ET[top] != NULL) edge[j][i].next = ET[top];
				ET[top] = &edge[j][i];

				polygon[i].e[j] = &edge[j][i];

			}

			if (ymax == -1 || ymin == HEIGHT) continue;
			polygon[i].dzx = - 2.0f * polygon[i].plane.x / (polygon[i].plane.z * WIDTH);
			polygon[i].dzy = 2.0f * polygon[i].plane.y / (polygon[i].plane.z * HEIGHT);
			polygon[i].dy = ymax - ymin;
			polygon[i].polygon_id = i;
			polygon[i].next = NULL;
			polygon[i].aet = NULL;
			polygon[i].apt = NULL;
			if (PT[ymax] != NULL) polygon[i].next = PT[ymax];
			PT[ymax] = &polygon[i];

		}

		active_polygon_table* APT_head = NULL;
		active_edge_table* AET_head = NULL;
		int count = 0;
		for (int i = HEIGHT - 1; i >= 0; i--) {

			for (int j = 0; j < WIDTH; j++) zbuffer[j] = 1000;

			polygon_table* current_PT = PT[i];
			while (current_PT != NULL) {

				active_polygon_table* new_APT = new active_polygon_table;
				new_APT->color = current_PT->color;
				new_APT->dy = current_PT->dy;
				new_APT->polygon_id = current_PT->polygon_id;
				new_APT->edges = 0;
				new_APT->next = NULL;
				if (APT_head != NULL) new_APT->next = APT_head;
				APT_head = new_APT;
				current_PT->apt = new_APT;

				current_PT = current_PT->next;

			}

			edge_table* current_ET = ET[i];
			while (current_ET != NULL) {

				if (!current_ET->out) {

					polygon_table* current_polygon = &polygon[current_ET->polygon_id];
					active_edge_table* current_AET;
					if (current_polygon->apt->edges == 0) {

						current_AET = new active_edge_table;
						current_AET->dxl = current_ET->dx;
						current_AET->dyl = current_ET->dy;
						current_AET->xl = current_ET->x;
						current_AET->polygon_id = current_ET->polygon_id;
						current_AET->zl = current_ET->z;
						current_AET->dzx = current_polygon->dzx;
						current_AET->dzy = current_polygon->dzy;
						current_AET->next = NULL;
						current_polygon->aet = current_AET;
						current_polygon->apt->edges++;
						if (AET_head != NULL) current_AET->next = AET_head;
						AET_head = current_AET;

					}
					else if (current_polygon->apt->edges == 1) {

						current_AET = current_polygon->aet;
						current_AET->dxr = current_ET->dx;
						current_AET->dyr = current_ET->dy;
						current_AET->xr = current_ET->x;
						current_polygon->apt->edges++;
						if (current_AET->xl > current_AET->xr || (current_AET->xl == current_AET->xr && current_AET->dxl > current_AET->dxr)) {

							swap(current_AET->xl, current_AET->xr);
							swap(current_AET->dxl, current_AET->dxr);
							swap(current_AET->dyl, current_AET->dyr);
							current_AET->zl = current_ET->z;

						}

					}
					else {

						current_AET = current_polygon->aet;
						if (current_AET->dyl == 0) {

							current_AET->dxl = current_ET->dx;
							current_AET->dyl = current_ET->dy;
							current_AET->xl = current_ET->x;
							current_AET->zl = current_ET->z;

						}
						else {

							current_AET->dxr = current_ET->dx;
							current_AET->dyr = current_ET->dy;
							current_AET->xr = current_ET->x;

						}

						if (current_AET->xl > current_AET->xr || (current_AET->xl == current_AET->xr && current_AET->dxl > current_AET->dxr)) {

							swap(current_AET->xl, current_AET->xr);
							swap(current_AET->dxl, current_AET->dxr);
							swap(current_AET->dyl, current_AET->dyr);

						}

					}
						
				}

				current_ET = current_ET->next;

			}

			active_polygon_table* current_APT = APT_head;
			active_polygon_table* last_APT = NULL;

			while (current_APT != NULL) {

				current_APT->dy--;
				if (current_APT->dy == 0) {

					if (current_APT == APT_head) {

						APT_head = current_APT->next;
						delete(current_APT);
						current_APT = APT_head;

					}
					else {

						last_APT->next = current_APT->next;
						delete(current_APT);
						current_APT = last_APT->next;

					}

				}
				else {

					last_APT = current_APT;
					current_APT = current_APT->next;

				}

			}

			active_edge_table* current_AET = AET_head;

			while (current_AET != NULL) {

				float z = current_AET->zl;
				int l = (int)current_AET->xl;
				int r = (int)current_AET->xr;

				if (l < 0) {

					z += -1.0f * l * current_AET->dzx;
					l = 0;

				}

				for (; l <= r && l <= WIDTH - 1; l++) {

					if (z < zbuffer[l] - EPSILON && z > -1.0f) {

						zbuffer[l] = z;
						PixelBuffer[(i * WIDTH + l) * 3] = polygon[current_AET->polygon_id].color.x;
						PixelBuffer[(i * WIDTH + l) * 3 + 1] = polygon[current_AET->polygon_id].color.y;
						PixelBuffer[(i * WIDTH + l) * 3 + 2] = polygon[current_AET->polygon_id].color.z;
						count++;

					}
					z += current_AET->dzx;
 
				}

				current_AET->dyl--;
				current_AET->dyr--;
				if (current_AET->dyl > 0) current_AET->xl += current_AET->dxl;
				if (current_AET->dyr > 0) current_AET->xr += current_AET->dxr;
				if (current_AET->dyl > 0) current_AET->zl += current_AET->dzx * current_AET->dxl + current_AET->dzy;
				current_AET = current_AET->next;

			}

			current_AET = AET_head;
			active_edge_table* last_AET = NULL;
			while (current_AET != NULL) {

				if (current_AET->dyl == 0 && current_AET->dyr == 0) {

					if (current_AET == AET_head) {

						AET_head = current_AET->next;
						delete(current_AET);
						current_AET = AET_head;

					}
					else {

						last_AET->next = current_AET->next;
						delete(current_AET);
						current_AET = last_AET->next;

					}

				}
				else {

					last_AET = current_AET;
					current_AET = current_AET->next;

				}

			}


		}

		if (count < 100) {

			bool x = true;

		}

		delete[](polygon);
		delete[](edge[0]);
		delete[](edge[1]);
		delete[](edge[2]);

		while (AET_head != NULL) {

			active_edge_table* current_AET = AET_head->next;
			delete(AET_head);
			AET_head = current_AET;

		}

		while (APT_head != NULL) {

			active_polygon_table* current_APT = APT_head->next;
			delete(APT_head);
			APT_head = current_APT;

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