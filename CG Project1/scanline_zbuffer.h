#pragma once
#ifndef SCANLINEZBUFFER
#define SCANLINEZBUFFER

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

struct active_polygon_table {

	int polygon_id;
	int dy;
	glm::vec3 color;
	int edges = 0;
	struct active_polygon_table* next;

};

struct active_edge_table {

	float xl, xr;
	float dxl, dxr;
	int dyl = 0, dyr = 0;
	float zl;
	float dzx;
	float dzy;
	int polygon_id;
	struct active_edge_table* next;
	struct active_polygon_table* active_polygon;
	bool AC_edge;

};

struct edge_table {

	float x;
	float dx;
	int dy;
	int polygon_id;
	struct edge_table* next;
	bool out = true;
	
};

struct polygon_table {

	vec4 plane;
	int polygon_id;
	int dy;
	float dzy;
	float dzx;
	float zl;
	glm::vec3 color;
	struct polygon_table* next; 
	struct edge_table* e[3];
	struct active_polygon_table* apt;
	struct active_edge_table* aet;

};

#endif // !SCANLINEZBUFFER
