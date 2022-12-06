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
#include "hierachical_zbuffer.h"

#include <omp.h>
#define NOMINMAX
#include <windows.h>  
# define	TIMING_BEGIN \
	{double tmp_timing_start = omp_get_wtime();

# define	TIMING_END(message) \
	{double tmp_timing_finish = omp_get_wtime();\
	double  tmp_timing_duration = tmp_timing_finish - tmp_timing_start;\
	printf("%s: %2.5f ms           \r", (message), tmp_timing_duration * 1000);}}


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
bool* lastVisibility;

unsigned char* PixelBuffer = new unsigned char[WIDTH * HEIGHT * 3];
vec3* color;

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
void hz_bottomRender(vec3 normal, vec3 top, vec3 downl, vec3 downr, hierachical_zbuffer* hz, unsigned char* pixelBuffer, vec3 color);
void hz_topRender(vec3 normal, vec3 top, vec3 downl, vec3 downr, hierachical_zbuffer* hz, unsigned char* pixelBuffer, vec3 color);
void ScanlineZbuffer();
void SimpleHierachicalZbuffer();
void HierachicalZbuffer();

int main() {

	load_obj_file("data/bunny_1k.obj");

	mesh.position = new vec3[position_queue.size()];
	mesh.ndc_position = new vec3[position_queue.size()];
	for (int i = 0; i < position_queue.size(); i++) mesh.position[i] = position_queue[i];

	mesh.normal = new vec3[normal_queue.size()];
	for (int i = 0; i < normal_queue.size(); i++) mesh.normal[i] = normal_queue[i];


	mesh.uv = new vec2[uv_queue.size()];
	for (int i = 0; i < uv_queue.size(); i++) mesh.uv[i] = uv_queue[i];

	mesh.triangle_indices = new triangle[indices_queue.size()];
	lastVisibility = new bool[indices_queue.size()];
	for (int i = 0; i < indices_queue.size(); i++) mesh.triangle_indices[i] = indices_queue[i], lastVisibility[i] = true;

	std::srand(std::time(0));
	color = new vec3[indices_queue.size()];
	for (int i = 0; i < indices_queue.size(); i++) {

		color[i].x = std::rand() % 255;
		color[i].y = std::rand() % 255;
		color[i].z = std::rand() % 255;

	}

	zbuffer = new float[WIDTH];

	//camera set up
	//cam = new camera(vec3(0.0f, 0.0f, -0.5f), vec3(0.0f, 1.0f, 0.0f), 90.0f, 25.0f);
	cam = new camera(vec3(-3.274809, -0.901280, -5.230831), vec3(0.0f, 1.0f, 0.0f), 55.899918, 24.500013);
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
	glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_HIDDEN);

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

		HierachicalZbuffer();
		//SimpleHierachicalZbuffer();
		//ScanlineZbuffer();

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

void ScanlineZbuffer() {

	polygon_table* polygon = new polygon_table[indices_queue.size()];
	edge_table* edge[3];
	edge[0] = new edge_table[indices_queue.size()];
	edge[1] = new edge_table[indices_queue.size()];
	edge[2] = new edge_table[indices_queue.size()];
	for (int i = 0; i < indices_queue.size(); i++) {

		int v[3];
		v[0] = mesh.triangle_indices[i].A();
		v[1] = mesh.triangle_indices[i].B();
		v[2] = mesh.triangle_indices[i].C();
		vec3* ndc_p = mesh.ndc_position;
		vec3 AB = vec3(ndc_p[v[0]] - ndc_p[v[1]]);
		vec3 CB = vec3(ndc_p[v[2]] - ndc_p[v[1]]);
		vec3 normal = cross(AB, CB);
		normal = normalize(normal);
		polygon[i].plane.x = normal.x;
		polygon[i].plane.y = normal.y;
		polygon[i].plane.z = normal.z;
		polygon[i].plane.w = -dot(normal, ndc_p[v[1]]);
		polygon[i].color.x = color[i].x;
		polygon[i].color.y = color[i].y;
		polygon[i].color.z = color[i].z;
		if (fabs(polygon[i].plane.z) < EPSILON) continue;

		int ymax = -1, ymin = HEIGHT;
		for (int j = 0; j < 3; j++) {

			edge[j][i].polygon_id = i;
			vec3 s = ndc_p[v[j]], t = ndc_p[v[(j + 1) % 3]];
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
		polygon[i].dzx = -2.0f * polygon[i].plane.x / (polygon[i].plane.z * WIDTH);
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

}

void SimpleHierachicalZbuffer() {

	int cutNum = 0;
	hierachical_zbuffer* hz = new hierachical_zbuffer(WIDTH, HEIGHT);
	for (int i = 0; i < indices_queue.size(); i++) {

		vec3 v[3];
		v[0] = mesh.ndc_position[mesh.triangle_indices[i].A()];
		v[1] = mesh.ndc_position[mesh.triangle_indices[i].B()];
		v[2] = mesh.ndc_position[mesh.triangle_indices[i].C()];
		vec3 AB = vec3(v[0] - v[1]);
		vec3 CB = vec3(v[2] - v[1]);
		vec3 normal = cross(AB, CB);
		normal = normalize(normal);
		if (fabs(normal.z) < EPSILON) continue;

		if (v[0].y < v[1].y) swap(v[0], v[1]);
		if (v[1].y < v[2].y) swap(v[1], v[2]);
		if (v[0].y < v[1].y) swap(v[0], v[1]);
		for (int i = 0; i < 3; i++) {

			v[i].x = (v[i].x + 1.0f) * WIDTH * 0.5f;
			v[i].y = (v[i].y + 1.0f) * HEIGHT * 0.5f;

		}

		vec2 minxy = vec2(v[0]), maxxy = vec2(v[0]);
		float nearest_distance = 1.0f;

		for (int i = 1; i < 3; i++) {

			if (v[i].x < minxy.x) minxy.x = v[i].x;
			if (v[i].y < minxy.y) minxy.y = v[i].y;
			if (v[i].x > maxxy.x) maxxy.x = v[i].x;
			if (v[i].y > maxxy.y) maxxy.y = v[i].y;
			if (v[i].z < nearest_distance) nearest_distance = v[i].z;

		}

		if (!hz->checkVisibilyty(minxy, maxxy, nearest_distance)) {

			cutNum++;
			continue;

		}

		if ((int)v[0].y == (int)v[1].y) hz_bottomRender(normal, v[0], v[1], v[2], hz, PixelBuffer, color[i]);
		else if ((int)v[1].y == (int)v[2].y) hz_topRender(normal, v[0], v[1], v[2], hz, PixelBuffer, color[i]);
		else {

			vec3 v3;
			float kx = -(v[0].x - v[2].x) / (v[0].y - v[2].y);
			float kz = -(v[0].z - v[2].z) / (v[0].y - v[2].y);
			v3.y = v[1].y;
			v3.x = v[0].x + kx * (v[0].y - v3.y);
			v3.z = v[0].z + kz * (v[0].y - v3.y);

			hz_bottomRender(normal, v[1], v3, v[2], hz, PixelBuffer, color[i]);
			hz_topRender(normal, v[0], v[1], v3, hz, PixelBuffer, color[i]);

		}

	}

	delete(hz);
	printf("Cut triangles: %d ", cutNum);

}

void HierachicalZbuffer() {

	vec3 AABB_v1 = vec3(1.0f), AABB_v2 = vec3(-1.0f);

	for (int i = 0; i < position_queue.size(); i++) {

		vec3 tmp = mesh.ndc_position[i];
		if (AABB_v1.x > tmp.x) AABB_v1.x = tmp.x;
		if (AABB_v1.y > tmp.y) AABB_v1.y = tmp.y;
		if (AABB_v1.z > tmp.z) AABB_v1.z = tmp.z;
		if (AABB_v2.x < tmp.x) AABB_v2.x = tmp.x;
		if (AABB_v2.y < tmp.y) AABB_v2.y = tmp.y;
		if (AABB_v2.z < tmp.z) AABB_v2.z = tmp.z;

	}

	if (AABB_v1.x < -1.0) AABB_v1.x = -1.0f;
	if (AABB_v1.y < -1.0) AABB_v1.y = -1.0f;
	if (AABB_v1.z < -1.0) AABB_v1.z = -1.0f;
	if (AABB_v2.x > 1.0) AABB_v2.x = 1.0f;
	if (AABB_v2.y > 1.0) AABB_v2.y = 1.0f;
	if (AABB_v2.z > 1.0) AABB_v2.z = 1.0f;

	int step = 1;
	int maxDeep = 0;
	while (step < WIDTH && step < HEIGHT) {

		step *= 2;
		maxDeep++;

	}

	int cutNum = 0;
	int tmp = 0;
	bool* currentVisibility = new bool[indices_queue.size()];
	memset(currentVisibility, 0, sizeof(bool) * indices_queue.size());
	octree* root = new octree(vec3(-1.0f), vec3(1.0f), maxDeep);
	for (int i = 0; i < indices_queue.size(); i++) 
		if (lastVisibility[i] == true) root->AddPolygon(i, mesh.ndc_position, mesh.triangle_indices);

	hierachical_zbuffer* hz = new hierachical_zbuffer(WIDTH, HEIGHT);
	root->Render(mesh, hz, cutNum, PixelBuffer, color, currentVisibility);
	root->Clear();
	delete(root);

	root = new octree(vec3(-1.0f), vec3(1.0f), maxDeep);
	for (int i = 0; i < indices_queue.size(); i++)
		if (lastVisibility[i] == false) root->AddPolygon(i, mesh.ndc_position, mesh.triangle_indices);

	//root->BuildTree(mesh.ndc_position, mesh.triangle_indices);
	root->Render(mesh, hz, cutNum, PixelBuffer, color, currentVisibility);
	root->Clear();
	delete(root);

	printf("Cut triangles: %d  ", cutNum);
	//memcpy(lastVisibility, currentVisibility, sizeof(bool) * indices_queue.size());
	
	for (int i = 0; i < indices_queue.size(); i++)
		lastVisibility[i] = currentVisibility[i];
	delete[](currentVisibility);
	delete(hz);

}