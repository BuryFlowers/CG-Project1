#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include "triangle.h"

using namespace glm;
using namespace std;

vector<vec3> position_queue;
vector<vec3> normal_queue;
vector<vec2> uv_queue;
vector<triangle> indices_queue;

void load_obj_file(const char* filename) {

	ifstream obj_file(filename);

	if (obj_file.fail()) {

		cerr << "Error: " << strerror(errno);
		exit(0);

	}
	const int MAX_LENGTH = 1024;
	char data[MAX_LENGTH] = { 0 };

	while (obj_file.getline(data, MAX_LENGTH)) {

		char type[MAX_LENGTH] = { 0 };
		if (sscanf(data, "%s", type) != 1) continue;

		if (strlen(type) == 1 && type[0] == 'v') {

			vec3 p;
			if (sscanf(data, "%s %f %f %f", type, &p.x, &p.y, &p.z) != 4) {

				printf("[Error] can't read potision from obj file!\n");
				exit(0);

			}

			position_queue.push_back(p);

		}

		else if (strlen(type) == 2 && type[0] == 'v' && type[1] == 'n') {

			vec3 n;
			if (sscanf(data, "%s %f %f %f", type, &n.x, &n.y, &n.z) != 4) {

				printf("[Error] can't read normal from obj file!\n");
				exit(0);

			}

			normal_queue.push_back(n);

		}

		else if (strlen(type) == 2 && type[0] == 'v' && type[1] == 't') {

			vec2 t;
			if (sscanf(data, "%s %f %f", type, &t.x, &t.y) != 3) {

				printf("[Error] can't read uvs from obj file!\n");
				exit(0);

			}

			uv_queue.push_back(t);

		}

		else if (strlen(type) == 1 && type[0] == 'f') {

			int p1, p2, p3;
			int n1, n2, n3;
			int t1, t2, t3;
			char c;
			if (sscanf(data, "%s %d %c %d %c %d %d %c %d %c %d %d %c %d %c %d", type, &p1, &c, &p2, &c, &p3, &n1, &c, &n2, &c, &n3, &t1, &c, &t2, &c, &t3) != 16) {

				printf("[Error] can't read indices from obj file!\n");
				exit(0);

			}

			indices_queue.push_back(triangle(p1, p2, p3, n1, n2, n3, t1, t2, t3));

		}

	}

	printf("[Success] an obj file has been loaded!\n");

}