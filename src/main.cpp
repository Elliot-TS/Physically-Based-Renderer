#include <iostream>
//#include <sstream> // For command-line arguments
#include "pbrt/samplers.h"
#include "pbrt/ray.h"
#include "pbrt/shapes.h"
#include "pbrt/primitive.h"
#include "pbrt/camera.h"
#include "pbrt/integrator.h"
#include "pbrt/util/display.h"

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
    // How to Read command line arguments
    //if (argc >= 2) {
        //std::istringstream iss( argv[1] );
        //if (!(iss >> val)) return -1;
    //} else return -1;

    // Ray tracing in a weekend
    int nx = 2000;
    int ny = 1000;
    int ns = 100;

    Uint32 c1 = Uint32(Color(100,100,100));
    std::cout << Color(c1) << std::endl;
    c1 = Color(c1) + Color(100,100,100);
    std::cout << Color(c1) << std::endl;
    c1 = (Color(c1) + Color(100,100,100));
    std::cout << Color(c1) << std::endl;
    c1 = (Color(c1) + Color(100,100,100));
    std::cout << Color(c1) << std::endl;

    Film film(nx, ny);
    Camera cam(&film);

    UniformSampler *sampler = new UniformSampler();

    Primitive *spheres[4] = { 
        new GeometricPrimitive(
                new Sphere(Point3f(0,0,-1), 0.5),
                new Lambertian(Vector3f(0.8,0.3,0.3), sampler)),
        new GeometricPrimitive(
                new Sphere(Point3f(0, -100.5, -1), 100),
                new Lambertian(Vector3f(0.8,0.8,0.01), sampler)),
        new GeometricPrimitive(
                new Sphere(Point3f(1,0,-1), 0.5),
                new Metal(Vector3f(0.8, 0.6, 0.2), 1.0, sampler)),
        new GeometricPrimitive(
                new Sphere(Point3f(-1,0,-1), 0.5),
                new Metal(Vector3f(0.8,0.8,0.8), 0.3, sampler))
    };
    SimpleAggregate sphPrims(spheres, 4);
    ImageTileIntegrator intr(&cam, sampler, &sphPrims, {});

    film.display->Open();
    intr.Render();
//    film.display->UpdateImage();


    // Convert command-line arguments to vector of strings TODO
    // Declare variables for parsed command line TODO
    // Process command-line arguments TODO
    // Initialize pbrt TODO
    // Parse provided scene description files TODO
    // Render the scene TODO
    // Clean up after rendering the scene TODO
}
