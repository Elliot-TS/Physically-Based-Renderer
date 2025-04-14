#include <iostream>
#include "pbrt/math/vecmath.h"
#include "pbrt/ray.h"

using namespace pbrt;
int main(int argc, char *argv[])
{
  Ray r(Point3f(0, 0, 0), Vector3f(1, 0, 0));
  std::cout << "Testing PBRT" << std::endl;
  std::cout << r << std::endl;
}
