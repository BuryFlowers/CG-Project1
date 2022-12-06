#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <iostream>
#include "camera.h"
#include "hierachical_zbuffer.h"

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

	//float ky = -(p1.x - p2.x) / (p1.y - p2.y);

}

bool shouldBeThrown(vec3 s, vec3 t) {

	if (s.y > MAX_HEIGHT && t.y > MAX_HEIGHT) return true;
	if (s.y < -1.0f && t.y < -1.0f) return true;
	if (s.x > -1.0f && s.x < MAX_WIDTH && s.y > -1.0f && s.y < MAX_HEIGHT) return false;
	if (t.x > -1.0f && t.x < MAX_WIDTH && t.y > -1.0f && t.y < MAX_HEIGHT) return false;
	if (t.z > -1.0f && t.z < 1.0f && t.z > -1.0f && t.z < 1.0f) return false;

	if (s.x > t.x) swap(s, t);

	vec3 direction = t - s;
	direction = normalize(direction);
	float len = length(t - s);

	float tnear, tfar;
	tnear = (-1.0f - s.x) / direction.x;
	tfar = (MAX_WIDTH - s.x) / direction.x;

	if (tnear > tfar) swap(tnear, tfar);

	float t1 = (-1.0f - s.y) / direction.y;
	float t2 = (MAX_HEIGHT - s.y) / direction.y;

	if (t1 > t2) swap(t1, t2);
	if (t1 > tnear) tnear = t1;
	if (t2 < tfar) tfar = t2;

	t1 = (-1.0f - s.z) / direction.z;
	t2 = (1.0f - s.z) / direction.z;

	if (t1 > t2) swap(t1, t2);
	if (t1 > tnear) tnear = t1;
	if (t2 < tfar) tfar = t2;

	if (tnear < tfar && tnear < len) return false;
	return false;

}

void hz_bottomRender(vec3 normal, vec3 topl, vec3 topr, vec3 down, hierachical_zbuffer* hz, unsigned char* pixelBuffer, vec3 color) {

	if (topl.x > topr.x) swap(topl, topr);

	if (topl.y > HEIGHT - 1 && down.y > HEIGHT - 1) return;
	if (topl.y < 0 && down.y < 0) return;
 
	float dxl = -(topl.x - down.x) / (topl.y - down.y);
	float dxr = -(topr.x - down.x) / (topr.y - down.y);
	float xl = topl.x;
	float xr = topr.x;
	float zl = topl.z;
	float dzx = -2.0f * normal.x / (normal.z * WIDTH);
	float dzy = 2.0f * normal.y / (normal.z * HEIGHT);
	float y = topl.y;
	int dy = (int)topl.y - (int)down.y;

	if  ((int)y > HEIGHT - 1) {

		int offset = (int)y - HEIGHT + 1;
		y = HEIGHT - 1;
		dy -= offset;
		xl += dxl * offset;
		xr += dxr * offset;
		zl += (dxl * dzx + dzy) * offset;

	}

	while(dy > 0 && y >= 0) {


		int l = (int)xl;
		int r = (int)xr;
		float z = zl;
		if (l < 0) z += -l * dzx, l = 0;
		r = std::min(WIDTH - 1, r);
		while (l <= r) {

			hz->UpdateDepth(l, y, z, pixelBuffer, color);
			z += dzx;
			l++;

		}

		dy--;
		y--;
		xl += dxl;
		xr += dxr;
		zl += dxl * dzx + dzy;


	}

}

void hz_topRender(vec3 normal, vec3 top, vec3 downl, vec3 downr, hierachical_zbuffer* hz, unsigned char* pixelBuffer,vec3 color) {

	if (downl.x > downr.x) swap(downl, downr);

	if (top.y > HEIGHT - 1 && downl.y > HEIGHT - 1) return;
	if (top.y < 0 && downl.y < 0) return;

	float dxl = -(top.x - downl.x) / (top.y - downl.y);
	float dxr = -(top.x - downr.x) / (top.y - downr.y);
	float xl = top.x;
	float xr = top.x;
	float zl = top.z;
	float dzx = -2.0f * normal.x / (normal.z * WIDTH);
	float dzy = 2.0f * normal.y / (normal.z * HEIGHT);
	float y = top.y;
	int dy = (int)top.y - (int)downl.y;

	if ((int)y > HEIGHT - 1) {

		int offset = (int)y - HEIGHT + 1;
		y = HEIGHT - 1;
		dy -= offset;
		xl += dxl * offset;
		xr += dxr * offset;
		zl += (dxl * dzx + dzy) * offset;

	}

	while (dy > 0 && y >= 0) {


		int l = (int)xl;
		int r = (int)xr;
		float z = zl;
		if (l < 0) z += -l * dzx, l = 0;
		r = std::min(WIDTH - 1, r);
		while (l <= r) {

			hz->UpdateDepth(l, y, z, pixelBuffer, color);
			z += dzx;
			l++;

		}

		dy--;
		y--;
		xl += dxl;
		xr += dxr;
		zl += dxl * dzx + dzy;


	}

}