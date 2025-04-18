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

bool LoadObject(
    std::string filename, Transform objectTransform,
    std::vector<Primitive *> *primitives, Sampler *sampler
)
{
  std::vector<std::unique_ptr<Shape>> triangles;

  objl::Loader Loader;
  bool loadout = Loader.LoadFile(filename);
  if (!loadout) return false;

  // Go through each loaded mesh
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
    TriangleMesh *mesh = new TriangleMesh(
        objectTransform, vertices, loadedMesh.Indices
    );
    Triangle::CreateTriangles(
        mesh, primitives,
        new Metal(
            Vector3f(
                loadedMesh.MeshMaterial.Kd.X,
                loadedMesh.MeshMaterial.Kd.Y,
                loadedMesh.MeshMaterial.Kd.Z
            ),
            0.8, sampler
        )
    );
  }
  return true;
}

// Based off RTW and PBRT
int main(int argc, char *argv[])
{
  // How to Read command line arguments
  std::string filename;
  if (argc >= 2) {
    filename = argv[1];
    std::cout << "Filename provided: " << filename << std::endl;
  }
  else {
    filename = "test.obj";
    std::cout << "No filename provided.  Using default: "
              << filename << std::endl;
  }


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
      &film, Point3f(3, 1, 2), Point3f(0, 0, 0), 50, aspectRatio
  );

  UniformSampler *sampler = new UniformSampler();

  std::vector<Primitive *> primitives = {
      new GeometricPrimitive(
          new Sphere(Point3f(0, -100.5, -1), 100),
          new Lambertian(Vector3f(0.34, 0.58, 0.47), sampler)
      ),
      new GeometricPrimitive(
          new Sphere(Point3f(-1.9, 0, 0), 0.1),
          new Metal(Vector3f(0.8, 0.6, 0.2), 1.0, sampler)
      ),
      new GeometricPrimitive(
          new Sphere(Point3f(-1, 0, -1), 0.5),
          new Metal(Vector3f(0.8, 0.8, 0.8), 0.3, sampler)
      ),
  };

  // Load objects
  bm.Start();
  LoadObject(
      filename, Rotate(45, Vector3f(0, 1, 0)), &primitives,
      sampler
  );
  std::cout << "Created Meshes: " << bm.GetTime() << "ms\n";

  bm.Start();
  BVHAggregate aggregate(
      primitives, primitives.size(),
      BVHAggregate::SplitMethod::SAH
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
