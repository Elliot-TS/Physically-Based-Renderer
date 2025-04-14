#include <iostream>
#include "pbrt/math/vecmath.h"
#include "pbrt/ray.h"

using namespace pbrt;
int main(int argc, char *argv[])
{
  std::cout << "Testing PBRT" << std::endl;

  Point3f origin(0, 0, -100);
  Point3f lookat(2.30, 0, -1);
  Ray r(origin, lookat - origin);
  Bounds3f b(Point3f(-1, -1, -1), Point3f(1, 1, 1));

  bool intersect = b.IntersectP(r.origin, r.direction);
  std::cout << (intersect ? "true" : "false") << std::endl;
}
