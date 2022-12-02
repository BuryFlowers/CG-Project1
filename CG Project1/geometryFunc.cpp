#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <iostream>
#include "camera.h"

using namespace glm;
using namespace std;

void attachToEdge(vec3 &p1, vec3 &p2) {

	float kx = -(p1.x - p2.x) / (p1.y - p2.y);
	float kz = -(p1.z - p2.z) / (p1.y - p2.y);
	float top = MAX_HEIGHT;
	if (p1.y > top) {

		p1.x = (p1.y - top) * kx + p1.x;
		p1.z = (p1.y - top) * kz + p1.z;
		p1.y = top;

	}
	else if (p1.y < -1.0f) {

		p1.x = (p2.y + 1.0f) * kx + p2.x;
		p1.z = (p2.y + 1.0f) * kz + p2.z;
		p1.y = -1.0f;

	}

	if (p2.y > top) {

		p2.x = (p2.y - top) * kx + p2.x;
		p2.z = (p2.y - top) * kz + p2.z;
		p2.y = top;

	}
	else if (p2.y < -1.0f) {

		p2.x = (p1.y + 1.0f) * kx + p1.x;
		p2.z = (p1.y + 1.0f) * kz + p1.z;
		p2.y = -1.0f;

	}

	//if (p1.y >= HEIGHT || p2.y >= HEIGHT) printf("\n%lf %lf\n", p1.y, p2.y);

}

bool shouldBeThrown(vec3 s, vec3 t) {


	if (s.x > -1.0f && s.x < MAX_WIDTH && s.y > -1.0f && s.y < MAX_HEIGHT) return false;
	if (t.x > -1.0f && t.x < MAX_WIDTH && t.y > -1.0f && t.y < MAX_HEIGHT) return false;

	if (s.x > t.x) swap(s, t);

	vec3 direction = t - s;
	direction = normalize(direction);
	float len = length(t - s);

	float tnear, tfar;
	tnear = (-1.0f - s.x) / direction.x;
	tfar = (MAX_WIDTH - s.x) / direction.x;

	float	t1 = (-1.0f - s.y) / direction.y;
	float	t2 = (MAX_HEIGHT - s.y) / direction.y;

	if (t1 > t2) swap(t1, t2);
	if (t1 > tnear) tnear = t1;
	if (t2 < tfar) tfar = t2;

	if (t1 < t2 && t1 < len) return false;
	return true;

}