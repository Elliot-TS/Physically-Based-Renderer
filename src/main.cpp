#include <iostream>
//#include <sstream> // For command-line arguments
#include "pbrt/samplers.h"
#include "pbrt/ray.h"
#include "pbrt/shapes.h"
#include "pbrt/primitive.h"
#include "pbrt/camera.h"
#include "pbrt/integrator.h"
#include "pbrt/util/display.h"
#include "pbrt/util/buffercache.h"

using namespace pbrt;

// Based off RTW

// Based off RTW
Float hit_sphere(const Point3f center, Float radius, const Ray& r) {
    Vector3f oc = r.origin - center;
    Float a = Dot(r.direction, r.direction);
    Float b = 2.0 * Dot(oc, r.direction);
    Float c = Dot(oc, oc) - Sqr(radius);
    Float discriminant = b*b - 4*a*c;
    if (discriminant < 0) return -1.0;
    else return (-b - sqrt(discriminant)) / (2.0*a);
}

// Based off RTW

// Based off RTW and PBRT
int main(int argc, char *argv[]) {
    using std::chrono::high_resolution_clock; using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;
    // How to Read command line arguments
    //if (argc >= 2) {
        //std::istringstream iss( argv[1] );
        //if (!(iss >> val)) return -1;
    //} else return -1;

    // Ray tracing in a weekend

    InitBufferCaches();
    Triangle::Init();

    float aspectRatio = 1.32;
    int ny = 1556;
    int nx = aspectRatio*ny;
    int ns = 100;

    Film film(nx, ny);
    Camera cam(&film, Point3f(0,1,-3), Vector3f(0,-0.6,1), 80, aspectRatio);

    UniformSampler *sampler = new UniformSampler();

    TriangleMesh mesh(
            {Point3f(-2,0,0), Point3f(0,0,0),
             Point3f(-2,2,0), Point3f(0,2,0)},
            {2,1,0, 1,2,3});
    auto triangles = Triangle::CreateTriangles(&mesh);

    const int numShapes = 6;
    Primitive *spheres[numShapes] = { 
        new GeometricPrimitive(
                new Sphere(Point3f(0,0,-1), 0.5),
                new Lambertian(Vector3f(0.8,0.3,0.3), sampler)),
        new GeometricPrimitive(
                new Sphere(Point3f(0, -100.5, -1), 100),
                new Lambertian(Vector3f(0.34,0.58,0.47), sampler)),
        new GeometricPrimitive(
                new Sphere(Point3f(1,0,-1), 0.5),
                new Metal(Vector3f(0.8, 0.6, 0.2), 1.0, sampler)),
        new GeometricPrimitive(
                new Sphere(Point3f(-1,0,-1), 0.5),
                new Metal(Vector3f(0.8,0.8,0.8), 0.3, sampler)),
        new GeometricPrimitive(
                triangles[0].get(),
                new Metal(Vector3f(0.8,0.8,0.8), 0.1, sampler)),
        new GeometricPrimitive(
                triangles[1].get(),
                new Metal(Vector3f(0.8,0.8,0.8), 0.1, sampler))
    };
    SimpleAggregate sphPrims(spheres, numShapes);
    ImageTileIntegrator intr(&cam, sampler, &sphPrims, {});

    film.display->Open();
    auto t1 = high_resolution_clock::now();
    intr.Render();
    auto t2 = high_resolution_clock::now();
    duration<double, std::milli> diff = t2 - t1;
    std::cout << "64 Samples" << std::endl;
    std::cout << "Time: " << diff.count()/1000  << " seconds" << std::endl;
    film.display->WaitUntilClosed();
//    film.display->UpdateImage();


    // Convert command-line arguments to vector of strings TODO
    // Declare variables for parsed command line TODO
    // Process command-line arguments TODO
    // Initialize pbrt TODO
    // Parse provided scene description files TODO
    // Render the scene TODO
    // Clean up after rendering the scene TODO
    /*
    */
}
