#version 330 core
layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 aColor;

out vec3 color;

uniform mat4 view;
uniform mat4 perspective;
//uniform vec3 polygonColor;

void main()
{
	gl_Position = perspective * view * vec4(aPos, 1.0);
	color = aColor;

}