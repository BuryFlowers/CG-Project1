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
#include "toojpeg.h"

#include <omp.h>
#define NOMINMAX
#include <windows.h>  
# define	TIMING_BEGIN \
	{double tmp_timing_start = omp_get_wtime();

# define	TIMING_END(message) \
	{double tmp_timing_finish = omp_get_wtime();\
	double  tmp_timing_duration = tmp_timing_finish - tmp_timing_start;\
	printf("%s: %2.5f ms           \n\n", (message), tmp_timing_duration * 1000);}}


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

triangle_mesh mesh;
polygon_table* PT[HEIGHT];
edge_table* ET[HEIGHT];
active_polygon_table* APT;
active_edge_table* AET;
float* zbuffer;
octree* root;

bool render = false;
unsigned char* PixelBuffer = new unsigned char[WIDTH * HEIGHT * 3];
vec3* color;
// output file
std::ofstream myFile;

float* vertices;
unsigned int* indices;

extern void load_obj_file(const char* filename);

//opengl callback
void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void mouse_callback(GLFWwindow* window, double posX, double posY);
void scroll_callback(GLFWwindow* window, double offsetX, double offsetY);
void processInput(GLFWwindow* window);

//functions
inline int getYCoord(float y) { return (int)((y + 1.0f) * 0.5f * HEIGHT + 0.5f); }
inline float getXCoord(float x) { return (x + 1.0f) * 0.5f * WIDTH; }
//cut edge coordinate to screen range
void attachToEdge(vec3& p1, vec3& p2);
//throw edges are out of screen
bool shouldBeThrown(vec3 s, vec3 t);
//split a triangle to two parts. One has a parallel bottom edge, the other one has a parallel top edge.
void hz_bottomRender(vec3 normal, vec3 top, vec3 downl, vec3 downr, hierachical_zbuffer* hz, unsigned char* pixelBuffer, vec3 color);
void hz_topRender(vec3 normal, vec3 top, vec3 downl, vec3 downr, hierachical_zbuffer* hz, unsigned char* pixelBuffer, vec3 color);
//three zbuffer alogorithms
void ScanlineZbuffer();
void SimpleHierachicalZbuffer();
void HierachicalZbuffer();

int main() {

	//load obj file
	load_obj_file("data/FullBody_Decimated_Small.obj");

	//merge the mesh data from seperate queues. Though normals and uvs are not used.
	mesh.position = new vec3[position_queue.size()];
	mesh.ndc_position = new vec3[position_queue.size()];
	for (int i = 0; i < position_queue.size(); i++) mesh.position[i] = position_queue[i];
	mesh.normal = new vec3[normal_queue.size()];
	for (int i = 0; i < normal_queue.size(); i++) mesh.normal[i] = normal_queue[i];
	mesh.uv = new vec2[uv_queue.size()];
	for (int i = 0; i < uv_queue.size(); i++) mesh.uv[i] = uv_queue[i];
	mesh.triangle_indices = new triangle[indices_queue.size()];
	for (int i = 0; i < indices_queue.size(); i++) mesh.triangle_indices[i] = indices_queue[i];

	//give every triangle a random color
	std::srand(std::time(0));
	color = new vec3[indices_queue.size()];
	for (int i = 0; i < indices_queue.size(); i++) {

		color[i].x = std::rand() % 255;
		color[i].y = std::rand() % 255;
		color[i].z = std::rand() % 255;

	}

	zbuffer = new float[WIDTH];

	//get the AABB of the mesh, then build the octree
	vec3 AABB_v1 = mesh.position[0], AABB_v2 = mesh.position[0];
	for (int i = 1; i < position_queue.size(); i++) {

		vec3 tmp = mesh.position[i];
		if (AABB_v1.x > tmp.x) AABB_v1.x = tmp.x;
		if (AABB_v1.y > tmp.y) AABB_v1.y = tmp.y;
		if (AABB_v1.z > tmp.z) AABB_v1.z = tmp.z;
		if (AABB_v2.x < tmp.x) AABB_v2.x = tmp.x;
		if (AABB_v2.y < tmp.y) AABB_v2.y = tmp.y;
		if (AABB_v2.z < tmp.z) AABB_v2.z = tmp.z;

	}

	//get octree's max deep/height
	int step = 1;
	int maxDeep = 0;
	while (step < WIDTH && step < HEIGHT) {

		step *= 2;
		maxDeep++;

	}
	root = new octree(AABB_v1, AABB_v2, maxDeep);
	for (int i = 0; i < indices_queue.size(); i++) root->AddPolygon(i, mesh.position, mesh.triangle_indices);
	printf("[Success] the octree has been built!\n\n");

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
	glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

	Shader shader("default.vs", "default.fs");

	//preapre vertices and indices to use opengl pipeline to render previewer 
	vertices = new float[position_queue.size() * 6];
	for (int i = 0; i < position_queue.size(); i++) {

		vertices[i * 6 + 0] = mesh.position[i].x;
		vertices[i * 6 + 1] = mesh.position[i].y;
		vertices[i * 6 + 2] = mesh.position[i].z;
		vertices[i * 6 + 3] = (rand() % 255) / 255.0f;
		vertices[i * 6 + 4] = (rand() % 255) / 255.0f;
		vertices[i * 6 + 5] = (rand() % 255) / 255.0f;

	}

	indices = new unsigned int[indices_queue.size() * 3];
	for (int i = 0; i < indices_queue.size(); i++) {

		indices[i * 3 + 0] = mesh.triangle_indices[i].A();
		indices[i * 3 + 1] = mesh.triangle_indices[i].B();
		indices[i * 3 + 2] = mesh.triangle_indices[i].C();

	}

	// init VBO VAO EBO
	unsigned int VBO, VAO, EBO;
	glGenVertexArrays(1, &VAO);
	glGenBuffers(1, &VBO);
	glGenBuffers(1, &EBO);

	glBindVertexArray(VAO);

	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(float) * position_queue.size() * 6, vertices, GL_STATIC_DRAW);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int) * indices_queue.size() * 3, indices, GL_STATIC_DRAW);

	// position attribute
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(0);
	// color attribute
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
	glEnableVertexAttribArray(1);

	//enable previewer's depth test
	glEnable(GL_DEPTH_TEST);

	while (!glfwWindowShouldClose(window)) {

		processInput(window);
		
		glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		//set coordinate space transform matrices
		mat4 viewMat = cam->GetViewMatrix();
		mat4 perspectiveMat = perspective(radians(cam->Zoom), RATIO, 0.1f, 100.0f);
		shader.setMat4("view", viewMat);
		shader.setMat4("perspective", perspectiveMat);

		//render the image if the user want
		if (render) {

			mat4 transform = perspectiveMat * viewMat;
			//get triangles' ndc
			for (int i = 0; i < position_queue.size(); i++) {

				vec4 tmp;
				tmp = transform * vec4(mesh.position[i], 1.0f);
				mesh.ndc_position[i] = vec3(tmp / fabs(tmp.w));

			}

			ScanlineZbuffer();
			SimpleHierachicalZbuffer();
			HierachicalZbuffer();
			render = false;

		}

		// render preview
		shader.use();
		glBindVertexArray(VAO);
		glDrawElements(GL_TRIANGLES, indices_queue.size() * 3, GL_UNSIGNED_INT, 0);
		glfwSwapBuffers(window);
		glfwPollEvents();

	}

	glfwTerminate();
	return 0;

}

// write a single byte compressed by tooJpeg
void myOutput(unsigned char byte)
{
	myFile << byte;
}
//use tooJpeg library to save image 
void SaveJPEG(const char* filename, unsigned char* PixelBuffer) {

	//flip the pixel buffer to get the right image
	for (int j = 0; j < HEIGHT / 2; j++) {

		int k = HEIGHT - 1 - j;
		for (int i = 0; i < WIDTH; i++) {
			
			swap(PixelBuffer[(k * WIDTH + i) * 3], PixelBuffer[(j * WIDTH + i) * 3]);
			swap(PixelBuffer[(k * WIDTH + i) * 3 + 1], PixelBuffer[(j * WIDTH + i) * 3 + 1]);
			swap(PixelBuffer[(k * WIDTH + i) * 3 + 2], PixelBuffer[(j * WIDTH + i) * 3 + 2]);

		}

	}

	//save 
	myFile = std::ofstream(filename, std::ios_base::out | std::ios_base::binary);
	const bool isRGB = true;  // true = RGB image, else false = grayscale
	const auto quality = 90;    // compression quality: 0 = worst, 100 = best, 80 to 90 are most often used
	const bool downsample = false; // false = save as YCbCr444 JPEG (better quality), true = YCbCr420 (smaller file)
	const char* comment = "result image"; // arbitrary JPEG comment
	auto ok = TooJpeg::writeJpeg(myOutput, PixelBuffer, WIDTH, HEIGHT, isRGB, quality, downsample, comment);

}

void ScanlineZbuffer() {

	//clear pixel buffer
	for (int i = 0; i < WIDTH * HEIGHT; i++) {

		PixelBuffer[i * 3] = (int)(255);
		PixelBuffer[i * 3 + 1] = (int)(255);
		PixelBuffer[i * 3 + 2] = (int)(255);

	}
	//clear polygon and edge table
	for (int i = 0; i < HEIGHT; i++) {

		PT[i] = NULL;
		ET[i] = NULL;

	}

	TIMING_BEGIN
	//allocate memory for every triangle and its edges
	polygon_table* polygon = new polygon_table[indices_queue.size()];
	edge_table* edge[3];
	edge[0] = new edge_table[indices_queue.size()];
	edge[1] = new edge_table[indices_queue.size()];
	edge[2] = new edge_table[indices_queue.size()];
	//traverse every triangle to create polygon and edge table
	for (int i = 0; i < indices_queue.size(); i++) {

		int v[3];
		v[0] = mesh.triangle_indices[i].A();
		v[1] = mesh.triangle_indices[i].B();
		v[2] = mesh.triangle_indices[i].C();
		vec3* ndc_p = mesh.ndc_position;
		//calculate the normal of the triangle's plane
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
		//if the plane is perpendicular to the xoy plane, it will not be rendered
		if (fabs(polygon[i].plane.z) < EPSILON) continue;

		//ymax is the max y of the triangle, and ymin is the min y of the triangle
		int ymax = -1, ymin = HEIGHT;
		//process three edges
		for (int j = 0; j < 3; j++) {

			edge[j][i].polygon_id = i;
			vec3 s = ndc_p[v[j]], t = ndc_p[v[(j + 1) % 3]];
			//make sure t is the higher vertex in y coordinate
			if (s.y > t.y) swap(s, t);
			//"out" means this edge is out of renderable edge range
			edge[j][i].out = false;
			if (shouldBeThrown(s, t)) {

				edge[j][i].out = true;
				continue;

			}
			//calculate the delta x when y decreases by 1
			edge[j][i].dx = -RATIO * (s.x - t.x) / (s.y - t.y);
			// use "attachToEdge" to restrain s and t's coordinate to valid screen space
			if (s.x < t.x) attachToEdge(s, t);
			else attachToEdge(t, s);

			//calculate the scan line number of the edge
			int top = getYCoord(t.y), down = getYCoord(s.y);
			edge[j][i].dy = top - down;
			edge[j][i].x = getXCoord(t.x);
			edge[j][i].z = t.z;
			//if the scan line number is 0, then there is no need to render it 
			if (edge[j][i].dy == 0) {

				edge[j][i].out = true;
				continue;

			}
			//update ymax and ymin
			if (top > ymax) ymax = top;
			if (down < ymin) ymin = down;

			//attach edge to edge table
			edge[j][i].next = NULL;
			if (ET[top] != NULL) edge[j][i].next = ET[top];
			ET[top] = &edge[j][i];

			//attach edge to polygon
			polygon[i].e[j] = &edge[j][i];

		}

		if (ymax == -1 || ymin == HEIGHT) continue;
		//dzx means delta z when x increases by 1
		polygon[i].dzx = -2.0f * polygon[i].plane.x / (polygon[i].plane.z * WIDTH);
		//dzy means delta z when y decreases by 1
		polygon[i].dzy = 2.0f * polygon[i].plane.y / (polygon[i].plane.z * HEIGHT);
		//calculate scan line number of this polygon
		polygon[i].dy = ymax - ymin;
		polygon[i].polygon_id = i;
		polygon[i].next = NULL;
		//clear its active edge table and active polygon table
		polygon[i].aet = NULL;
		polygon[i].apt = NULL;
		//attach this polygon to polygon table
		if (PT[ymax] != NULL) polygon[i].next = PT[ymax];
		PT[ymax] = &polygon[i];

	}

	active_polygon_table* APT_head = NULL;
	active_edge_table* AET_head = NULL;
	//start scan line rendering
	for (int i = HEIGHT - 1; i >= 0; i--) {

		//clear zbuffer
		for (int j = 0; j < WIDTH; j++) zbuffer[j] = 1000;

		//get every polygon in this scan line
		polygon_table* current_PT = PT[i];
		while (current_PT != NULL) {

			//attach every polygon to a new active polygon table
			active_polygon_table* new_APT = new active_polygon_table;
			new_APT->color = current_PT->color;
			new_APT->dy = current_PT->dy;
			new_APT->polygon_id = current_PT->polygon_id;

			//set number of edge that has been attached to a active edge table 
			new_APT->edges = 0;

			//attach the new active polygon to active polyon table
			new_APT->next = NULL;
			if (APT_head != NULL) new_APT->next = APT_head;
			APT_head = new_APT;
			//set this polygon's active polygon 
			current_PT->apt = new_APT;

			current_PT = current_PT->next;

		}

		//get every edge in this scan line
		edge_table* current_ET = ET[i];
		while (current_ET != NULL) {

			//if current edge may need to render
			if (!current_ET->out) {

				polygon_table* current_polygon = &polygon[current_ET->polygon_id];
				active_edge_table* current_AET;
				//if the polgon of this edge has no active edge
				if (current_polygon->apt->edges == 0) {

					//create a new active edge
					current_AET = new active_edge_table;
					//attach the edge to the left side of the new active edge
					current_AET->dxl = current_ET->dx;
					current_AET->dyl = current_ET->dy;
					current_AET->xl = current_ET->x;
					current_AET->polygon_id = current_ET->polygon_id;
					current_AET->zl = current_ET->z;
					current_AET->dzx = current_polygon->dzx;
					current_AET->dzy = current_polygon->dzy;
					current_AET->next = NULL;
					//set polygon's active edge
					current_polygon->aet = current_AET;
					//the polygon's active edge increases by 1
					current_polygon->apt->edges++;
					//attach the new active edge to active edge table
					if (AET_head != NULL) current_AET->next = AET_head;
					AET_head = current_AET;

				}
				//if the polgon of this edge already has one active edge
				else if (current_polygon->apt->edges == 1) {

					//attach the edge to the right side of the its active edge
					current_AET = current_polygon->aet;
					current_AET->dxr = current_ET->dx;
					current_AET->dyr = current_ET->dy;
					current_AET->xr = current_ET->x;
					current_polygon->apt->edges++;
					//swap the left side and the right side if needed
					if (current_AET->xl > current_AET->xr || (current_AET->xl == current_AET->xr && current_AET->dxl > current_AET->dxr)) {

						swap(current_AET->xl, current_AET->xr);
						swap(current_AET->dxl, current_AET->dxr);
						swap(current_AET->dyl, current_AET->dyr);
						current_AET->zl = current_ET->z;

					}

				}
				//if the polgon of this edge already has two active edge
				else {

					current_AET = current_polygon->aet;
					//find the edge should be thrown
					if (current_AET->dyl == 0) {

						//attach the edge to the left side of the its active edge
						current_AET->dxl = current_ET->dx;
						current_AET->dyl = current_ET->dy;
						current_AET->xl = current_ET->x;
						current_AET->zl = current_ET->z;

					}
					else {

						//attach the edge to the right side of the its active edge
						current_AET->dxr = current_ET->dx;
						current_AET->dyr = current_ET->dy;
						current_AET->xr = current_ET->x;

					}

					//swap the left side and the right side if needed
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
		//throw edges should be thrown
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
		//start render every active edge
		while (current_AET != NULL) {

			float z = current_AET->zl;
			int l = (int)current_AET->xl;
			int r = (int)current_AET->xr;
			//if l is out of left side of screen, correct the l and z
			if (l < 0) {

				z += -1.0f * l * current_AET->dzx;
				l = 0;

			}
			//start to render
			for (; l <= r && l <= WIDTH - 1; l++) {

				//if thie edge is the nearer one
				if (z < zbuffer[l] - EPSILON && z > -1.0f) {

					zbuffer[l] = z;
					PixelBuffer[(i * WIDTH + l) * 3] = polygon[current_AET->polygon_id].color.x;
					PixelBuffer[(i * WIDTH + l) * 3 + 1] = polygon[current_AET->polygon_id].color.y;
					PixelBuffer[(i * WIDTH + l) * 3 + 2] = polygon[current_AET->polygon_id].color.z;

				}
				//change the z
				z += current_AET->dzx;

			}

			//decrease dy by 1 after finishing rendering this scan line
			current_AET->dyl--;
			current_AET->dyr--;
			//correct x and z
			if (current_AET->dyl > 0) current_AET->xl += current_AET->dxl;
			if (current_AET->dyr > 0) current_AET->xr += current_AET->dxr;
			if (current_AET->dyl > 0) current_AET->zl += current_AET->dzx * current_AET->dxl + current_AET->dzy;
			//next active edge
			current_AET = current_AET->next;

		}

		current_AET = AET_head;
		active_edge_table* last_AET = NULL;
		//delete active polygon that has been fully rendered
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
	
	//release memory
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

	TIMING_END("Scan line time")
	SaveJPEG("./result/scanline_zbuffer.jpeg", PixelBuffer);
	
}

void SimpleHierachicalZbuffer() {

	//clear pixel buffer
	for (int i = 0; i < WIDTH * HEIGHT; i++) {

		PixelBuffer[i * 3] = (int)(255);
		PixelBuffer[i * 3 + 1] = (int)(255);
		PixelBuffer[i * 3 + 2] = (int)(255);

	}

	TIMING_BEGIN
	//use "cutNum" to record the number of triangles that are rejected by hierachical zbuffer
	int cutNum = 0;
	hierachical_zbuffer* hz = new hierachical_zbuffer(WIDTH, HEIGHT);
	//traverse every triangle
	for (int i = 0; i < indices_queue.size(); i++) {

		vec3 v[3];
		v[0] = mesh.ndc_position[mesh.triangle_indices[i].A()];
		v[1] = mesh.ndc_position[mesh.triangle_indices[i].B()];
		v[2] = mesh.ndc_position[mesh.triangle_indices[i].C()];
		//calculate the normal of this triangle's plane
		vec3 AB = vec3(v[0] - v[1]);
		vec3 CB = vec3(v[2] - v[1]);
		vec3 normal = cross(AB, CB);
		normal = normalize(normal);
		//if the plane is perpendicular to the xoy plane, it will not be rendered
		if (fabs(normal.z) < EPSILON) continue;

		//make v[0] be the highest vertex on y coordinate
		if (v[0].y < v[1].y) swap(v[0], v[1]);
		if (v[1].y < v[2].y) swap(v[1], v[2]);
		if (v[0].y < v[1].y) swap(v[0], v[1]);
		for (int i = 0; i < 3; i++) {

			//map ndc to window space 
			v[i].x = (v[i].x + 1.0f) * WIDTH * 0.5f;
			v[i].y = (v[i].y + 1.0f) * HEIGHT * 0.5f;

		}

		//get AABB and its nearest distance towards camera
		vec2 minxy = vec2(v[0]), maxxy = vec2(v[0]);
		float nearest_distance = 1.0f;
		for (int i = 1; i < 3; i++) {

			if (v[i].x < minxy.x) minxy.x = v[i].x;
			if (v[i].y < minxy.y) minxy.y = v[i].y;
			if (v[i].x > maxxy.x) maxxy.x = v[i].x;
			if (v[i].y > maxxy.y) maxxy.y = v[i].y;
			if (v[i].z < nearest_distance) nearest_distance = v[i].z;

		}

		//check if it has potential to be rendered
		if (!hz->checkVisibilyty(minxy, maxxy, nearest_distance)) {

			cutNum++;
			continue;

		}

		//render the triangle
		//if this triangle doesn't need to be split to two parts
		if ((int)v[0].y == (int)v[1].y) hz_bottomRender(normal, v[0], v[1], v[2], hz, PixelBuffer, color[i]);
		else if ((int)v[1].y == (int)v[2].y) hz_topRender(normal, v[0], v[1], v[2], hz, PixelBuffer, color[i]);
		//if this triangle needs to be split to two parts
		else {

			//calculate the middle edge's another vertex
			vec3 v3;
			float kx = -(v[0].x - v[2].x) / (v[0].y - v[2].y);
			float kz = -(v[0].z - v[2].z) / (v[0].y - v[2].y);
			v3.y = v[1].y;
			v3.x = v[0].x + kx * (v[0].y - v3.y);
			v3.z = v[0].z + kz * (v[0].y - v3.y);
			//reder two parts
			hz_bottomRender(normal, v[1], v3, v[2], hz, PixelBuffer, color[i]);
			hz_topRender(normal, v[0], v[1], v3, hz, PixelBuffer, color[i]);

		}

	}

	delete(hz);
	
	printf("Cut triangles: %d \n", cutNum);
	TIMING_END("Simple hierachical zbuffer time")
	SaveJPEG("./result/simple_hierachical_zbuffer.jpeg", PixelBuffer);

}

void HierachicalZbuffer() {

	//clear pixel buffer
	for (int i = 0; i < WIDTH * HEIGHT; i++) {

		PixelBuffer[i * 3] = (int)(255);
		PixelBuffer[i * 3 + 1] = (int)(255);
		PixelBuffer[i * 3 + 2] = (int)(255);

	}

	mat4 transform = perspective(radians(cam->Zoom), RATIO, 0.1f, 100.0f) * cam->GetViewMatrix();
	TIMING_BEGIN

	//use "cutNum" to record the number of triangles that are rejected by hierachical zbuffer
	int cutNum = 0;
	hierachical_zbuffer* hz = new hierachical_zbuffer(WIDTH, HEIGHT);
	//render the triangles that are not rejected by hierachical zbuffer
	root->Render(mesh, hz, cutNum, PixelBuffer, color, transform, true);
	//render the triangles that are rejected by hierachical zbuffer
	root->Render(mesh, hz, cutNum, PixelBuffer, color, transform, false);
	delete(hz);

	printf("Cut triangles: %d \n", cutNum);
	TIMING_END("Hierachical zbuffer time")

	SaveJPEG("./result/hierachical_zbuffer.jpeg", PixelBuffer);

}