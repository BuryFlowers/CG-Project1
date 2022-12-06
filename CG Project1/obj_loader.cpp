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

			int p[3] = { -1 };
			int n[3] = { -1 };
			int t[3] = { -1 };
			char s[3][MAX_LENGTH];
			if (sscanf(data, "%s %s %s %s", type, s[0], s[1], s[2]) != 4) {

				printf("[Error] can't read indices from obj file!\n");
				exit(0);

			}

			else for (int i = 0; i < 3; i++) {

					int tmp = 0, j = 0;
					while (s[i][j] != '/' && s[i][j] != ' ' && j < strlen(s[i])) {

						tmp *= 10;
						tmp += s[i][j] - '0';
						j++;

					}

					p[i] = tmp - 1;
					j++;
					tmp = 0;

					while (s[i][j] != '/' && s[i][j] != ' ' && j < strlen(s[i])) {

						tmp *= 10;
						tmp += s[i][j] - '0';
						j++;

					}

					if (s[i][j] == '/' && j < strlen(s[i])) n[i] = tmp - 1, j++;
					tmp = 0;

					while (s[i][j] != '/' && s[i][j] != ' ' && j < strlen(s[i])) {

						tmp *= 10;
						tmp += s[i][j] - '0';
						j++;

					}

					if (s[i][j] == '/' && j < strlen(s[i])) t[i] = tmp - 1, j++;
					tmp = 0;

			}

			indices_queue.push_back(triangle(p[0], p[1], p[2], n[0], n[1], n[2], t[0], t[1], t[2]));

		}

	}

	printf("[Success] an obj file has been loaded!\n\n");

}