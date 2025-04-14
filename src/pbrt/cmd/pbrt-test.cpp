#include <iostream>
#include "pbrt/math/vecmath.h"
#include "pbrt/ray.h"
#include "pbrt/shapes.h"

using namespace pbrt;

void test_sphere_intersection(
    Point3f origin, Point3f lookat, Sphere s,
    bool expectsIntersectBounds, bool expectsIntersectSphere
)
{
  std::cout << "Testing sphere:\n\tcenter: " << s.center
            << "\n\tradius: " << s.radius << std::endl
            << "Intersection with ray at\n\torigin: " << origin
            << "\n\tlookat: " << lookat << std::endl;

  Ray r(origin, Normalize(lookat - origin));

  Bounds3f b = s.Bounds();

  std::cout << (b.IntersectP(r.origin, r.direction)
                    ? "Ray intersects bounds"
                    : "Ray does not intersect bounds")
            << " Expected = " << expectsIntersectBounds
            << std::endl;
  std::cout << (s.Intersect(r)
                    ? "Ray intersects Sphere"
                    : "Ray does not intersect Sphere")
            << " Expected = " << expectsIntersectSphere
            << std::endl;
}

int main(int argc, char *argv[])
{
  std::cout << "Testing PBRT" << std::endl;
  test_sphere_intersection(
      Point3f(0, 0, 100), Point3f(1.9, 0, 0),
      Sphere(Point3f(0, 0, 0), 1), false, false
  );

  test_sphere_intersection(
      Point3f(0, 0, 100), Point3f(-1.9, 0, 0),
      Sphere(Point3f(0, 0, 0), 1), false, false
  );

  test_sphere_intersection(
      Point3f(0, 0, 100), Point3f(-0.99, 0, 0),
      Sphere(Point3f(0, 0, 0), 1), true, true
  );

  test_sphere_intersection(
      Point3f(0, 0, 100), Point3f(1.0, 0, 0),
      Sphere(Point3f(0, 0, 0), 1), true, false
  );
}
