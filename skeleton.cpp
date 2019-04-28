#include "skeleton.h"
#include "SDLauxiliary.h"
#include "model_triangle.h"
#include <SDL.h>
#include <glm/glm.hpp>
#include <iostream>
#include <utility>

using namespace std;
using glm::vec3;
using glm::ivec2;
using glm::mat3;

const int SCREEN_WIDTH = 500;
const int SCREEN_HEIGHT = 500;
SDL_Surface *screen;
int t; // time
vector<Triangle> triangles;

// camera parameters
const float focalLength = 500;
vec3 cameraPos(0, 0, -3.00);

// light parameters
vec3 lightPos(0, -0.5, -0.7);
vec3 lightPower = 10.1f * vec3(1, 1, 1);
vec3 indirectLightPowerPerArea = 0.5f * vec3(1, 1, 1);

// rotation variables
mat3 R(1.0f);
float yaw = 0;

vec3 rightt;
vec3 down;
vec3 forward;

vec3 currentColor;

vec3 currentNormal;
vec3 currentReflectance;

float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH];

static void PixelShader(const Pixel &p) {
  int x = p.x;
  int y = p.y;
  if (p.zinv > depthBuffer[y][x]) {

    depthBuffer[y][x] = p.zinv;

    vec3 dir(lightPos.x - p.pos3d.x,
             lightPos.y - p.pos3d.y,
             lightPos.z - p.pos3d.z);

	float matrix_dot_product = glm::dot(glm::normalize(dir),
										glm::normalize(currentNormal));

	float max_dot_product = std::max(matrix_dot_product, 0.0f);

	// Sphere 4*dir^2*pi
	float surface_area = static_cast<float>(4 * pow(glm::length(dir), 2) * M_PI);

    vec3 D = lightPower * (max_dot_product / surface_area);

    vec3 color =
        currentReflectance * (D + indirectLightPowerPerArea) * currentColor;

    PutPixelSDL(screen, x, y, color);
  }
}

static void DrawPolygon(const vector<Vertex> &vertices)
{
  int V = vertices.size();
  vector<Pixel> vertexPixels(V);

  for (int i = 0; i < V; ++i)
    VertexShader(vertices[i], vertexPixels[i]);

  vector<Pixel> leftPixels(0);
  vector<Pixel> rightPixels(0);

  ComputePolygonRows(vertexPixels, leftPixels, rightPixels);
  DrawPolygonRows(leftPixels, rightPixels);
}


int main(int argc, char *argv[]) {
  updateMatrix();
  // Get the triangles to be drawn
  LoadTestModel(triangles);
  screen = InitializeSDL(SCREEN_WIDTH, SCREEN_HEIGHT);
  // Set start value for timer.
  t = SDL_GetTicks();

  while (NoQuitMessageSDL()) {
    Update();
    Draw();
  }

  SDL_SaveBMP(screen, "screenshot.bmp");
  return 0;
}

void Update() {
  // Compute frame time:
  int t2 = SDL_GetTicks();
  float dt = float(t2 - t);
  t = t2;

  cout << "Render time: " << dt << " ms." << endl;

  Uint8 *keystate = SDL_GetKeyState(0);
  if (keystate[SDLK_UP]) {
    cameraPos += forward * 0.2f;
  } else if (keystate[SDLK_DOWN]) {
    cameraPos -= forward * 0.2f;
  } else if (keystate[SDLK_RIGHT]) {
    yaw += 0.01;
    updateMatrix();
  } else if (keystate[SDLK_LEFT]) {
    yaw -= 0.01;
    updateMatrix();
  } else if (keystate[SDLK_RSHIFT]) {
    cameraPos += 0.4f * down;
  } else if (keystate[SDLK_RCTRL]) {
    cameraPos -= 0.4f * down;
  } else if (keystate[SDLK_z]) {
    lightPos += 0.4f * forward;
  } else if (keystate[SDLK_s]) {
    lightPos -= 0.4f * forward;
  } else if (keystate[SDLK_d]) {
    lightPos += 0.4f * rightt;
  } else if (keystate[SDLK_q]) {
    lightPos -= 0.4f * rightt;
  } else if (keystate[SDLK_r]) {
    lightPos -= 0.4f * down;
  } else if (keystate[SDLK_f]) {
    lightPos += 0.4f * down;
  }

}

/*
 * Recompute the matrices for moving the camera
 * as they might be updated by the user
 * Update the rotation matrix
 */
void updateMatrix() {

  // Rotation matrice depending on yaw
  vec3 c1(cos(yaw), 0, -sin(yaw));
  vec3 c2(0       , 1, 0        );
  vec3 c3(sin(yaw), 0, cos(yaw) );

  /*
   *  cos  0  -sin
   *  0    1   0
   *  sin  0   cos
   */
  R = mat3(c1, c2, c3);
  rightt  = vec3(R[0][0], R[0][1], R[0][2]);
  down    = vec3(R[1][0], R[1][1], R[1][2]);
  forward = vec3(R[2][0], R[2][1], R[2][2]);
}

void Draw() {
  for (int y = 0; y < SCREEN_HEIGHT; ++y) {
    for (int x = 0; x < SCREEN_WIDTH; ++x) {
      depthBuffer[y][x] = 0;
    }
  }

  // Create rectangle of screen size which will be the futur screen
  SDL_FillRect(screen, 0, 0);


  if (SDL_MUSTLOCK(screen)) {
    SDL_LockSurface(screen);
  }

  // For all triangles
  for (int i = 0; i < triangles.size(); ++i) {
    // Get its color, normalization vector and position
    currentColor = triangles[i].color;
    currentNormal = triangles[i].normal;
    currentReflectance = vec3(1, 1, 1);

    vector<Vertex> vertices(3);
    vertices[0].position = triangles[i].v0;
    vertices[1].position = triangles[i].v1;
    vertices[2].position = triangles[i].v2;

    DrawPolygon(vertices);
    // DrawPolygonEdges(vertices);
  }

  if (SDL_MUSTLOCK(screen)) {
    SDL_UnlockSurface(screen);
  }

  // Swap with new screen
  SDL_UpdateRect(screen, 0, 0, 0, 0);
}

void VertexShader(const Vertex &v, Pixel &p)
{
  // Compute the coordinate of the point taking care of the position and angle
  // of the camera
  vec3 P = (v.position - cameraPos) * R;

  p.x = static_cast<int>(focalLength * P.x / P.z) + SCREEN_WIDTH / 2;
  p.y = static_cast<int>(focalLength * P.y / P.z) + SCREEN_HEIGHT / 2;
  p.zinv = 1 / P.z;
  p.pos3d = v.position;

  // vec3
  // dir(lightPos.x-v.position.x,lightPos.y-v.position.y,lightPos.z-v.position.z);

  // vec3 D = lightPower *
  // (std::max(glm::dot(glm::normalize(dir),glm::normalize(v.normal)),static_cast<float>(0))/static_cast<float>(4*pow(glm::length(dir),2)*M_PI));

  // p.illumination = v.reflectance *
  // (D+indirectLightPowerPerArea)*currentColor;
}

void Interpolate(Pixel a, Pixel b, vector<Pixel> &result)
{
  // Get the number of pixels to be drawn
  int N = result.size();
  // Find the step between the two points depending on the number of points
  // the line will draw
  //        0      1       2       3    
  // 0                          B(0,3)
  // 1                      ----|
  // 2                  ----|
  // 3              ----|
  // 4          ----|
  // 5      ----|
  // 6   A(6,0)
  //
  // N = max(|A.x - B.x|, |A.y - B.y|) = max(6,3) = 6
  // step = (B - A) / (N-1) = (-6/6,3/6) = (-1,1/2)
  //   So decrement x by 1 and increment y by 1/2

  vec3 step = vec3(vec3(b.x, b.y, b.zinv) - vec3(a.x, a.y, a.zinv)) /
                   float(max(N - 1, 1));
  // vec3 stepIllum = vec3(b.illumination - a.illumination) /
  // float(max(N-1,1));
  vec3 stepPos3d =
      vec3(b.pos3d * b.zinv - a.pos3d * a.zinv) / float(max(N - 1, 1));

  vec3 current(a.x, a.y, a.zinv);
  // vec3 curentIllum(a.illumination);
  vec3 curentPos3d(a.pos3d * a.zinv);
  for (int i = 0; i < N; ++i) {
    result[i].x = current.x;
    result[i].y = current.y;
    result[i].zinv = current.z;
    // result[i].illumination = curentIllum;
    result[i].pos3d = curentPos3d / current.z;

    current += step;
    curentPos3d += stepPos3d;
    // curentIllum += stepIllum;
  }
}

void DrawPolygonEdges(const vector<Vertex> &vertices) {
  int V = vertices.size();
  // Transform each vertex from 3D world position to 2D
  // image position:
  vector<Pixel> projectedVertices(V);
  for (int i = 0; i < V; ++i) {
    VertexShader(vertices[i], projectedVertices[i]);
  }

  // Loop over all vertices and draw the edge from it
  // to the next vertex:
  for (int i = 0; i < V; ++i) {
    int j = (i + 1) % V; // The next vertex
    vec3 color(1, 1, 1);

    vector<Pixel> line = lineFromInterpolation(projectedVertices[i],
                                               projectedVertices[j]);

    for (int x = 0; x < line.size(); ++x) {
        if (line[x].x >= 0 &&
            line[x].x < SCREEN_WIDTH &&
            line[x].y >= 0 &&
            line[x].y < SCREEN_HEIGHT) {

            PutPixelSDL(screen, line[x].x, line[x].y, color);
        }
    }
  }
}

void ComputePolygonRows(const vector<Pixel> &vertexPixels,
                        vector<Pixel> &leftPixels,
                        vector<Pixel> &rightPixels)
{

  // Find max and min y-value of the polygon and compute the number of rows it
  // occupies.
  int minY = +numeric_limits<int>::max();
  int maxY = -numeric_limits<int>::max();

  for (int i = 0; i < vertexPixels.size(); ++i) {
    if (vertexPixels[i].y > maxY) {
      maxY = vertexPixels[i].y;
    }

    if (vertexPixels[i].y < minY) {
      minY = vertexPixels[i].y;
    }
  }

  int occupied_rows = maxY - minY + 1;

  leftPixels.resize(occupied_rows);
  rightPixels.resize(occupied_rows);

  for (int i = 0; i < occupied_rows; ++i) {
    leftPixels[i].x  = numeric_limits<int>::max();
    rightPixels[i].x = numeric_limits<int>::min();

    leftPixels[i].y  = minY + i;
    rightPixels[i].y = minY + i;
  }

  // Loop through all edges of the polygon and use linear interpolation to find
  // the x-coordinate for each row it occupies. Update the corresponding
  // values in rightPixels and leftPixels.
  int size = vertexPixels.size();
  for (int i = 0; i < size; ++i) {
    vector<Pixel> edge =
        lineFromInterpolation(vertexPixels[i], vertexPixels[(i + 1) % size]);

    for (int j = 0; j < edge.size(); ++j) {
      int row = edge[j].y - minY;
      if (row >= 0 && row < leftPixels.size()) {
        if (edge[j].x < leftPixels[row].x) {
          leftPixels[row] = edge[j];
        }
        if (edge[j].x > rightPixels[edge[j].y - minY].x) {
          rightPixels[row] = edge[j];
        }
      }
    }
  }
}

vector<Pixel> lineFromInterpolation(Pixel a, Pixel b) {
  // Compute |a - b|
  ivec2 delta = glm::abs(ivec2(a.x, a.y) - ivec2(b.x, b.y));
  // Find the max number of pixels that separate the two points
  int pixels = glm::max(delta.x, delta.y) + 1;
  // Create a line of the max number of pixel found
  vector<Pixel> line(pixels);
  // Initialize the line
  Interpolate(a, b, line);

  return line;
}

void DrawPolygonRows(const vector<Pixel> &leftPixels,
                     const vector<Pixel> &rightPixels)
{

  for (int y = 0; y < leftPixels.size(); ++y) {
    vector<Pixel> line = lineFromInterpolation(leftPixels[y], rightPixels[y]);

    for (int x = 0; x < line.size(); ++x) {
      if (line[x].x >= 0 &&
            line[x].x < SCREEN_WIDTH &&
            line[x].y >= 0 &&
          line[x].y < SCREEN_HEIGHT) {

          PixelShader(line[x]);
        }
    }
  }
}
