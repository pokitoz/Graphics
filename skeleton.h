#ifndef SKELETON_H_
#define SKELETON_H_

#include "SDLauxiliary.h"
#include "model_triangle.h"
#include <glm/glm.hpp>


struct Pixel {
  int x;
  int y;
  float zinv;
  glm::vec3 pos3d;
};

struct Vertex {
  glm::vec3 position;
};

void Update();
void Draw();

void DrawPolygonEdges(const std::vector<Vertex> &vertices);
void updateMatrix();

void Interpolate(Pixel a, Pixel b, std::vector<Pixel> &result);

void ComputePolygonRows(const std::vector<Pixel> &vertexPixels,
                        std::vector<Pixel> &leftPixels,
                        std::vector<Pixel> &rightPixels);

std::vector<Pixel> lineFromInterpolation(Pixel a, Pixel b);

void DrawPolygonRows(const std::vector<Pixel> &leftPixels,
                     const std::vector<Pixel> &rightPixels);

void VertexShader(const Vertex &v, Pixel &p);

#endif /* SKELETON_H_ */
