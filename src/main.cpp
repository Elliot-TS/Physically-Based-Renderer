#include <nanobench.h>
#include <iostream>
#include "ObjLoader/OBJ_Loader.h"
#include "pbrt/aggregates.h"
#include "pbrt/camera.h"
#include "pbrt/integrator.h"
#include "pbrt/primitive.h"
#include "pbrt/ray.h"
#include "pbrt/samplers.h"
#include "pbrt/shapes.h"
#include "pbrt/util/benchmark.h"
#include "pbrt/util/buffercache.h"
#include "pbrt/util/display.h"

using namespace pbrt;

// Based off RTW

// Based off RTW
Float hit_sphere(
    const Point3f center, Float radius, const Ray &r
)
{
  Vector3f oc = r.origin - center;
  Float a = Dot(r.direction, r.direction);
  Float b = 2.0 * Dot(oc, r.direction);
  Float c = Dot(oc, oc) - Sqr(radius);
  Float discriminant = b * b - 4 * a * c;
  if (discriminant < 0)
    return -1.0;
  else
    return (-b - sqrt(discriminant)) / (2.0 * a);
}

// Based off RTW

// Based off RTW and PBRT
int main(int argc, char *argv[])
{
  // How to Read command line arguments
  // if (argc >= 2) {
  // std::istringstream iss( argv[1] );
  // if (!(iss >> val)) return -1;
  //} else return -1;

  // Ray tracing in a weekend

  Benchmark bm;

  InitBufferCaches();
  Triangle::Init();

  float aspectRatio = 1.32;
  int ny = 1556;
  int nx = aspectRatio * ny;
  int ns = 100;

  Film film(nx, ny);
  Camera cam(
      &film, Point3f(0, 1, -3), Vector3f(0, -0.6, 1), 80,
      aspectRatio
  );

  UniformSampler *sampler = new UniformSampler();

  bm.Start();
  objl::Loader Loader;
  bool loadout = Loader.LoadFile("test2.obj");
  if (!loadout) return 1;

  // Go through each loaded mesh
  std::vector<TriangleMesh> triMeshes;
  std::vector<std::unique_ptr<Shape>> triangles;
  for (int i = 0; i < Loader.LoadedMeshes.size(); ++i) {
    objl::Mesh loadedMesh = Loader.LoadedMeshes[i];

    // Vertices
    std::vector<Point3f> vertices;
    for (int j = 0; j < loadedMesh.Vertices.size(); ++j) {
      vertices.push_back(Point3f(
          loadedMesh.Vertices[j].Position.X,
          loadedMesh.Vertices[j].Position.Y,
          loadedMesh.Vertices[j].Position.Z
      ));
    }

    // Triangle Mesh
    triMeshes.push_back(
        TriangleMesh(vertices, loadedMesh.Indices)
    );
    auto newTriangles =
        Triangle::CreateTriangles(&triMeshes.back());
    std::move(
        newTriangles.begin(), newTriangles.end(),
        std::back_inserter(triangles)
    );
  }

  std::vector<Primitive *> primitives = {
      new GeometricPrimitive(
          new Sphere(Point3f(0, -100.5, -1), 100),
          new Lambertian(Vector3f(0.34, 0.58, 0.47), sampler)
      ),
      new GeometricPrimitive(
          new Sphere(Point3f(1, 0, -1), 0.5),
          new Metal(Vector3f(0.8, 0.6, 0.2), 1.0, sampler)
      ),
      new GeometricPrimitive(
          new Sphere(Point3f(-1, 0, -1), 0.5),
          new Metal(Vector3f(0.8, 0.8, 0.8), 0.3, sampler)
      ),
  };
  int numShapes = 3;
  for (int i = 0; i < triangles.size(); ++i) {
    primitives.push_back(new GeometricPrimitive(
        triangles[i].get(),
        new Metal(Vector3f(0.8, 0.1, 0.2), 0.3, sampler)
    ));
    numShapes++;
  }
  std::cout << "Created Meshes: " << bm.GetTime() << " ms \n";

  bm.Start();
  BVHAggregate aggregate(
      primitives, numShapes, BVHAggregate::SplitMethod::SAH
  );
  // SimpleAggregate aggregate(&primitives[0], numShapes);
  std::cout << "Create aggregate: " << bm.GetTime() << " ms\n";

  ImageTileIntegrator intr(&cam, sampler, &aggregate, {});

  film.display->Open();

  bm.Start();
  intr.Render();
  std::cout << "64 Samples" << std::endl;
  std::cout << "Time: " << bm.GetTime(true) << " seconds"
            << std::endl;
  film.display->WaitUntilClosed();
  //    film.display->UpdateImage();


  // Convert command-line arguments to vector of strings TODO
  // Declare variables for parsed command line TODO
  // Process command-line arguments TODO
  // Initialize pbrt TODO
  // Parse provided scene description files TODO
  // Render the scene TODO
  // Clean up after rendering the scene TODO
}
