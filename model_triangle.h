#ifndef MODEL_TRIANGLE_H
#define MODEL_TRIANGLE_H

#include <glm/glm.hpp>
#include <vector>
#include <math.h>

using glm::vec3;

class Triangle
{
public:
	vec3 v0;
	vec3 v1;
	vec3 v2;
	vec3 normal;
	vec3 color;

	Triangle(vec3 v0, vec3 v1, vec3 v2, vec3 color)
		: v0(v0), v1(v1), v2(v2), color(color)
	{
		ComputeNormal();
	}

	void ComputeNormal()
	{
		vec3 e1 = v1-v0;
		vec3 e2 = v2-v0;
		normal = glm::normalize(glm::cross(e2, e1));
	}
};

vec3 red(    0.80f, 0.10f, 0.10f );
vec3 white(  0.80f, 0.80f, 0.80f );

float remap_custom(float value,
				   float low1,
				   float high1,
				   float low2,
				   float high2)
{
	if(value > high1) {
		value = high1;
	}

	return low2 + (value - low1) * (high2 - low2) / (high1 - low1);
}

void DrawFace(std::vector<Triangle>& triangles,
			  vec3 p1, vec3 p2, vec3 p3, vec3 p4,
			  vec3 color) {

	triangles.push_back(Triangle(p1,p2,p3,color));
	triangles.push_back(Triangle(p4,p2,p1,color));

}

void DrawCircle(std::vector<Triangle>& triangles,
				vec3 center, float diameter, vec3 color)
{
	float r = diameter / 2.0f;

	int total_longitude = 100;
	int total_latitude  = 100;
	
	vec3 flat_globe[total_longitude+1][total_latitude + 1];
	
	
	for(int i = 0; i < total_longitude + 1; i++) {
	
		float longitude = remap_custom(i,
							           0, total_longitude,
									   -M_PI/2, M_PI/2);

		for(int j = 0; j < total_latitude + 1; j++) {

			float latitude = remap_custom(j,
										  0, total_latitude,
										  -M_PI, M_PI);

			float x = r * sin(longitude) * cos(latitude);
			float y = r * sin(longitude) * sin(latitude);
			float z = r * cos(longitude);

			vec3 a(center.x + x,
				   center.y + y,
				   center.z + z);

			flat_globe[i][j] = a;
		}
	}

	for(int i = 0; i < total_longitude; i++) {
		for(int j = 0; j < total_latitude; j++) {

				// A B
				// C D
				vec3 A(flat_globe[i][j]);
				vec3 B(flat_globe[i][j+1]);
				vec3 C(flat_globe[i+1][j]);
				vec3 D(flat_globe[i+1][j+1]);

				triangles.push_back(Triangle(A, D, B, white));
				triangles.push_back(Triangle(A, D, C, white));
		}
	}
}



void DrawBox_coordinate(std::vector<Triangle>& triangles,
						vec3 A, vec3 B, vec3 C, vec3 D,
						vec3 E, vec3 F, vec3 G, vec3 H)
{
	/* 
	 * The cube looks like:
	 *	  G----H
	 *	 /    /|
	 *	E----F |
	 *	| |  | |
	 *	| C----D
	 *	|/   |/
	 *	A----B
     */

	// Face (ABEF)
	// Triangle EBA and EBF
	DrawFace(triangles, E,B,A,F, white);

	// Face (BDHF)
	// Triangle FDH and FDB
	DrawFace(triangles, F,D,H,B, white);

	// Face (CDHG)
	// Triangle HCD and HCG
	DrawFace(triangles, H,C,D,G, white);

	// Face (ACGE)
	// Triangle ECA and ECG
	DrawFace(triangles, E,C,A,G, white);

	// Face (EFHG)
	// Triangle GFE and GFH
	DrawFace(triangles, G,F,E,H, white);	

	// Face (ABCD)
	// Triangle ABC and ABD
	DrawFace(triangles, A,B,C,D, white);
}

void DrawBox_size(std::vector<Triangle>& triangles,
			 	  vec3 base,
			 	  float size_x, float size_y, float size_z)
{

	float x_left = base.x;
	float y_down = base.y;		
	float z_back = base.z;

	float x_right = x_left + size_x;
	float y_up = y_down + size_y;
	float z_front = z_back + size_z;

	vec3 A(x_left , y_down, z_front);
	vec3 B(x_right, y_down, z_front);
	vec3 C(x_left , y_down, z_back);
	vec3 D(x_right, y_down, z_back);

	vec3 E(x_left , y_up  , z_front);
	vec3 F(x_right, y_up  , z_front);
	vec3 G(x_left , y_up  , z_back);
	vec3 H(x_right, y_up  , z_back);

	DrawBox_coordinate(triangles, A,B,C,D,E,F,G,H);
}

void LoadTestModel(std::vector<Triangle>& triangles)
{
	triangles.clear();
	triangles.reserve(20*20*5);

	float L = 555;

	vec3 A(L,0,0);
	vec3 B(0,0,0);
	vec3 C(L,0,L);
	vec3 D(0,0,L);

	vec3 E(L,L,0);
	vec3 F(0,L,0);
	vec3 G(L,L,L);
	vec3 H(0,L,L);

	// Floor:
	triangles.push_back(Triangle(C, B, A, white));
	triangles.push_back(Triangle(C, D, B, white));

	// Left wall
	triangles.push_back(Triangle(A, E, C, white));
	triangles.push_back(Triangle(C, E, G, white));

	// Right wall
	triangles.push_back(Triangle(F, B, D, white));
	triangles.push_back(Triangle(H, F, D, white));

	// Ceiling
	triangles.push_back(Triangle(E, F, G, white));
	triangles.push_back(Triangle(F, H, G, white));

	vec3 center(255,255,255);
	DrawCircle(triangles, center, 100, white);

	// Scale to the volume [-1,1]^3
	for( size_t i=0; i<triangles.size(); ++i )
	{
		triangles[i].v0 *= 2/L;
		triangles[i].v1 *= 2/L;
		triangles[i].v2 *= 2/L;

		triangles[i].v0 -= vec3(1,1,1);
		triangles[i].v1 -= vec3(1,1,1);
		triangles[i].v2 -= vec3(1,1,1);

		triangles[i].v0.x *= -1;
		triangles[i].v1.x *= -1;
		triangles[i].v2.x *= -1;

		triangles[i].v0.y *= -1;
		triangles[i].v1.y *= -1;
		triangles[i].v2.y *= -1;

		triangles[i].ComputeNormal();
	}
}

#endif
