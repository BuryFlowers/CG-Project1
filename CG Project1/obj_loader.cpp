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

	//read the file
	ifstream obj_file(filename);

	if (obj_file.fail()) {

		cerr << "Error: " << strerror(errno);
		exit(0);

	}
	//set max length of every line
	const int MAX_LENGTH = 1024;
	char data[MAX_LENGTH] = { 0 };

	while (obj_file.getline(data, MAX_LENGTH)) {

		char type[MAX_LENGTH] = { 0 };
		//get the line's first word
		if (sscanf(data, "%s", type) != 1) continue;

		//if the word is 'v', this line will give the pisition of a vertex
		if (strlen(type) == 1 && type[0] == 'v') {

			vec3 p;
			//read the position
			if (sscanf(data, "%s %f %f %f", type, &p.x, &p.y, &p.z) != 4) {

				printf("[Error] can't read potision from obj file!\n");
				exit(0);

			}

			//push it to queue
			position_queue.push_back(p);

		}

		//if the word is 'vn', this line will give the normal of a vertex
		else if (strlen(type) == 2 && type[0] == 'v' && type[1] == 'n') {

			vec3 n;
			//read the normal
			if (sscanf(data, "%s %f %f %f", type, &n.x, &n.y, &n.z) != 4) {

				printf("[Error] can't read normal from obj file!\n");
				exit(0);

			}

			//push it to queue
			normal_queue.push_back(n);

		}

		//if the word is 'vt', this line will give the texture coordinatess of a vertex
		else if (strlen(type) == 2 && type[0] == 'v' && type[1] == 't') {

			vec2 t;
			//read the uv
			if (sscanf(data, "%s %f %f", type, &t.x, &t.y) != 3) {

				printf("[Error] can't read uvs from obj file!\n");
				exit(0);

			}

			//push it to queue
			uv_queue.push_back(t);

		}

		//if the word is 'f', this line will give the indices of a triangle
		else if (strlen(type) == 1 && type[0] == 'f') {

			int p[3] = { 0 };
			int n[3] = { 0 };
			int t[3] = { 0 };
			char s[3][MAX_LENGTH];
			//read three strings
			if (sscanf(data, "%s %s %s %s", type, s[0], s[1], s[2]) != 4) {

				printf("[Error] can't read indices from obj file!\n");
				exit(0);

			}

			else for (int i = 0; i < 3; i++) {
					
					//get the position index
					int tmp = 0, j = 0;
					while (s[i][j] != '/' && s[i][j] != ' ' && j < strlen(s[i])) {

						tmp *= 10;
						tmp += s[i][j] - '0';
						j++;

					}

					p[i] = tmp - 1;

					//get the normal index
					j++;
					tmp = 0;
					while (s[i][j] != '/' && s[i][j] != ' ' && j < strlen(s[i])) {

						tmp *= 10;
						tmp += s[i][j] - '0';
						j++;

					}

					if (s[i][j] == '/' && j < strlen(s[i])) n[i] = tmp - 1, j++;
					
					//get the uv index
					tmp = 0;
					while (s[i][j] != '/' && s[i][j] != ' ' && j < strlen(s[i])) {

						tmp *= 10;
						tmp += s[i][j] - '0';
						j++;

					}

					if (s[i][j] == '/' && j < strlen(s[i])) t[i] = tmp - 1, j++;
					tmp = 0;

			}

			//push it to the queue
			indices_queue.push_back(triangle(p[0], p[1], p[2], n[0], n[1], n[2], t[0], t[1], t[2]));

		}

	}

	printf("[Success] an obj file has been loaded!\n");

}