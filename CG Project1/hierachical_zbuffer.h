#pragma once
#ifndef HIERACHICAL_ZBUFFER
#define HIERACHICAL_ZBUFFER

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <vector>
#include <algorithm>
#include "triangle.h"

using namespace glm;

class hierachical_zbuffer {

public:

	hierachical_zbuffer(int width, int height) {
		
		int length = 1;
		max_level = 0;
		while (length < width && length < height) {

			length *= 2;
			max_level++;

		}
		max_level--;

		level_width = new int[max_level + 1];
		level_height = new int[max_level + 1];
		length = 1;
		for (int i = 0; i <= max_level; i++) {

			level_width[i] = width / length + 1;
			level_height[i] = height / length + 1;
			zbuffer[i] = new float[level_width[i] * level_height[i]];
			renderCount[i] = new int[level_width[i] * level_height[i]];
			length *= 2;

		}

		for (int i = 0; i <= max_level; i++)
			for (int j = 0; j < level_width[i]; j++)
				for (int k = 0; k < level_height[i]; k++) {

					zbuffer[i][k * level_width[i] + j] = 2.0f;
					renderCount[i][k * level_width[i] + j] = 0;

				}

		w = width;
		h = height;

	}
	~hierachical_zbuffer() {

		delete[](level_width);
		delete[](level_height);
		for (int i = 0; i <= max_level; i++) delete[](zbuffer[i]);
		for (int i = 0; i <= max_level; i++) delete[](renderCount[i]);

	}

	int MaxLevel() { return max_level; }
	void UpdateDepth(int x, int y, float depth, unsigned char* pixelBuffer, vec3 color) {

		if (depth > -1.0f && depth < zbuffer[0][y * w + x]) {

			bool firstObject = false;
			if (zbuffer[0][y * w + x] == 2.0f) firstObject = true;
			zbuffer[0][y * w + x] = depth;
			pixelBuffer[(y * w + x) * 3] = color.x;
			pixelBuffer[(y * w + x) * 3 + 1] = color.y;
			pixelBuffer[(y * w + x) * 3 + 2] = color.z;

			
			int step = 2;
			for (int i = 1; i <= max_level; i++) {

				int level_x = x / step;
				int level_y = y / step;
				int index = level_x + level_y * level_width[i];

				if (depth > zbuffer[i][index] || zbuffer[i][index] == 2.0f) zbuffer[i][index] = depth;

				if (firstObject == true) renderCount[i][index]++;
				step *= 2;

			}

		}

	}
	bool checkVisibilyty(vec2 minxy, vec2 maxxy, float nearest_distance) {

		if (minxy.x >= w || minxy.y >= h) return false;
		else if (maxxy.x < 0 || maxxy.y < 0) return false;

		if (minxy.x < 0 || maxxy.x >= w) return true;
		else if (minxy.y < 0 || maxxy.y >= h) return true;

		int step = 1;
		for (int i = 1; i <= max_level; i++) step *= 2;
		for (int i = max_level; i >= 1; i--) {

			int lx, ly, rx, ry;
			lx = (int)minxy.x / step;
			ly = (int)minxy.y / step;
			rx = (int)maxxy.x / step;
			ry = (int)maxxy.y / step;

			//if (lx + 1 < rx || rx + 1 < ry) break;

			int fullBlockNum = 0;
			for (int x = lx; x <= rx; x++)
				for (int y = ly; y <= ry; y++) {

					int block_width = step, block_height = step;

					if (x == level_width[i] - 1) block_width = (w - 1) % step + 1;
					if (y == level_height[i] - 1) block_height = (h - 1) % step + 1;
					int index = x + y * level_width[i];

					if (renderCount[i][index] == block_width * block_height) {

						if (zbuffer[i][index] > nearest_distance) return true;
						fullBlockNum++;

					}

				}

			if (fullBlockNum == (rx - lx + 1) * (ry - ly + 1)) return false;

			step /= 2;

		}

		return true;

	}

private:
	int w;
	int h;
	int max_level;
	int* level_width;
	int* level_height;
	float* zbuffer[20];
	int* renderCount[20];

};

void hz_bottomRender(vec3 normal, vec3 topl, vec3 topr, vec3 down, hierachical_zbuffer* hz, unsigned char* pixelBuffer, vec3 color);
void hz_topRender(vec3 normal, vec3 top, vec3 downl, vec3 downr, hierachical_zbuffer* hz, unsigned char* pixelBuffer, vec3 color);

class octree {

public:

	octree(vec3 v1, vec3 v2, int max_deep) {

		for (int i = 0; i < 8; i++) children[i] = NULL;
		polygons.clear();
		visibility.clear();
		center = (v1 + v2) * 0.5f;
		length = center - v1;
		maxDeep = max_deep;

	}

	void AABB(mat4& transform, vec3& minxy, vec3& maxxy) {

		for (int i = 0; i < 8; i++) {

			vec4 v = vec4(center, 1.0f);
			v.x += direction_x[i] * length.x;
			v.y += direction_y[i] * length.y;
			v.z += direction_z[i] * length.z;
			v = transform * v;
			v /= fabs(v.w);
				
			v.x = (v.x + 1.0f) * WIDTH * 0.5f;
			v.y = (v.y + 1.0f) * HEIGHT * 0.5f;
			if (minxy.x > v.x || i == 0) minxy.x = v.x;
			if (minxy.y > v.y || i == 0) minxy.y = v.y;
			if (minxy.z > v.z || i == 0) minxy.z = v.z;
			if (maxxy.x < v.x || i == 0) maxxy.x = v.x;
			if (maxxy.y < v.y || i == 0) maxxy.y = v.y;
			if (maxxy.z < v.z || i == 0) maxxy.z = v.z;

		}

	}

	void AddPolygon(int x, vec3* position, triangle* triangles) {

		polygonNum++;
		if (maxDeep == 0) {

			polygons.push_back(x);
			visibility.push_back(true);
			return;

		}
		vec3 v[3];
		v[0] = position[triangles[x].A()];
		v[1] = position[triangles[x].B()];
		v[2] = position[triangles[x].C()];

		vec3 AABBv1 = v[0], AABBv2 = v[0];
		for (int j = 1; j < 3; j++) {

			if (AABBv1.x > v[j].x) AABBv1.x = v[j].x;
			if (AABBv1.y > v[j].y) AABBv1.y = v[j].y;
			if (AABBv1.z > v[j].z) AABBv1.z = v[j].z;
			if (AABBv2.x < v[j].x) AABBv2.x = v[j].x;
			if (AABBv2.y < v[j].y) AABBv2.y = v[j].y;
			if (AABBv2.z < v[j].z) AABBv2.z = v[j].z;

		}

		for (int j = 0; j < 8; j++) {

			vec3 box_v1 = center;
			vec3 box_v2 = center;
			box_v2.x += length.x * direction_x[j];
			box_v2.y += length.y * direction_y[j];
			box_v2.z += length.z * direction_z[j];
			if (box_v1.x > box_v2.x) std::swap(box_v1.x, box_v2.x);
			if (box_v1.y > box_v2.y) std::swap(box_v1.y, box_v2.y);
			if (box_v1.z > box_v2.z) std::swap(box_v1.z, box_v2.z);

			if (box_v1.x <= AABBv1.x && AABBv1.x <= box_v2.x && box_v1.y <= AABBv1.y && AABBv1.y <= box_v2.y && box_v1.z <= AABBv1.z && AABBv1.z <= box_v2.z) {

				if (box_v1.x <= AABBv2.x && AABBv2.x <= box_v2.x && box_v1.y <= AABBv2.y && AABBv2.y <= box_v2.y && box_v1.z <= AABBv2.z && AABBv2.z <= box_v2.z) {

					if (children[j] == NULL) children[j] = new octree(box_v1, box_v2, maxDeep - 1);
					children[j]->AddPolygon(x, position, triangles);
					return;

				}
				

			}

		}

		polygons.push_back(x);
		visibility.push_back(true);

	}

	void Render(triangles& mesh, hierachical_zbuffer* hz, int& cutNum, unsigned char* PixelBuffer, vec3* color, mat4& transform, bool mode) {

		vec3 minxy, maxxy;
		this->AABB(transform, minxy, maxxy);
	
		bool flag = hz->checkVisibilyty(vec2(minxy), vec2(maxxy), center.z - length.z);

		if (!flag) {

			cutNum += polygonNum;
			return;

		}

		for (int i = 0; i < 8; i++)
			if (children[i] != NULL)
				for (int j = i + 1; j < 8; j++) 
					if (children[j] != NULL) {

						vec3 iminxy, imaxxy, jminxy, jmaxxy;
						children[i]->AABB(transform, iminxy, imaxxy);
						children[j]->AABB(transform, jminxy, jmaxxy);
						if (iminxy.z > jminxy.z) std::swap(children[i], children[j]);

					}

		for (int i = 0; i < 4; i++)
			if (children[i] != NULL) children[i]->Render(mesh, hz, cutNum, PixelBuffer, color, transform, mode);


		for (int x = 0; x < polygons.size(); x++) {

			if (visibility[x] != mode) continue;
			int i = polygons[x];
			visibility[x] = false;
			vec3 v[3];
			v[0] = mesh.ndc_position[mesh.triangle_indices[i].A()];
			v[1] = mesh.ndc_position[mesh.triangle_indices[i].B()];
			v[2] = mesh.ndc_position[mesh.triangle_indices[i].C()];
			vec3 AB = vec3(v[0] - v[1]);
			vec3 CB = vec3(v[2] - v[1]);
			vec3 normal = cross(AB, CB);
			normal = normalize(normal);
			if (fabs(normal.z) < EPSILON) return;

			if (v[0].y < v[1].y) std::swap(v[0], v[1]);
			if (v[1].y < v[2].y) std::swap(v[1], v[2]);
			if (v[0].y < v[1].y) std::swap(v[0], v[1]);
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

			visibility[x] = true;

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

		for (int i = 4; i < 8; i++)
				if (children[i] != NULL) children[i]->Render(mesh, hz, cutNum, PixelBuffer, color, transform, mode);

	}

	void Clear() {

		for (int i = 0; i < 8; i++)
			if (children[i] != NULL) {

				children[i]->Clear();
				delete(children[i]);

			}
		
		polygons.clear();
		visibility.clear();

	}

private:
	int direction_x[8] = {-1, 1, -1, 1, -1, 1, -1, 1};
	int direction_y[8] = {-1, -1, 1, 1, -1, -1, 1, 1};
	int direction_z[8] = {-1, -1, -1, -1, 1, 1, 1, 1};
	octree* children[8];
	vec3 center;
	vec3 length;
	int maxDeep;
	int polygonNum = 0;
	std::vector<int> polygons;
	std::vector<bool> visibility;

};
#endif // !HIERACHICAL_ZBUFFER
