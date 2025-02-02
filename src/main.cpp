#include <iostream>
//#include <sstream> // For command-line arguments
#include "pbrt/samplers.h"
#include "pbrt/ray.h"
#include "pbrt/shapes.h"
#include "pbrt/cpu/primitive.h"
#include "pbrt/camera.h"

using namespace pbrt;


// Based off RTW
Point3f random_in_unit_sphere(UniformSampler &s)
{
    Vector3f p;
    do {
        p = 2.0*s.s_Vector3f() - Vector3f(1,1,1);
    } while (LengthSquared(p) >= 1.0);

    return Point3f(p);
}

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
Vector3f color(const Ray& r, const Primitive& shape, int bounces) {
    auto si = shape.Intersect(r);
    if (si) {
        Ray scattered;
        Vector3f attenuation;
        if (bounces > 0 && 
                si->interaction.material->scatter(
                    r, si->interaction, attenuation, scattered)
           )
        {
            return HorizontalProduct(
                    attenuation,
                    color(scattered, shape, bounces+1));
        }
        else return Vector3f(0,0,0);
    }

    Vector3f unit_direction= Normalize(r.direction);
    Float t = 0.5*(unit_direction.y + 1.0);
    return (1.0-t)  * Vector3f(1,1,1) 
            + t     * Vector3f(0.5, 0.7, 1.0);
}

// Based off RTW and PBRT
int main(int argc, char *argv[]) {
    // How to Read command line arguments
    //if (argc >= 2) {
        //std::istringstream iss( argv[1] );
        //if (!(iss >> val)) return -1;
    //} else return -1;

    // Ray tracing in a weekend
    int nx = 200;
    int ny = 100;
    int ns = 100;
    std::cout << "P3\n" << nx << " " << ny << "\n255\n";

    Camera cam;

    UniformSampler *sampler = new UniformSampler();

    Primitive *spheres[2] = { 
        new GeometricPrimitive(
                new Sphere(Point3f(0,0,-1), 0.5),
                new Lambertian(Vector3f(0.8,0.3,0.3), sampler)),
        new GeometricPrimitive(
                new Sphere(Point3f(0, -100.5, -1), 100),
                new Lambertian(Vector3f(0.8,0.8,0.8), sampler)),
    };
    SimpleAggregate sphPrims(spheres, 2);


    for (int j = ny-1; j >= 0; j--) {
        for (int i = 0; i < nx; i++) {
            Vector3f col(0,0,0);

            for (int s = 0; s < ns; s++) {
                Float u = Float(i + sampler->sample()) / Float(nx);
                Float v = Float(j + sampler->sample()) / Float(ny);

                Ray r = cam.get_ray(u, v);
                col += color(r, sphPrims, 4);
            }
            col /= ns;
            col = Vector3f(
                    std::sqrt(col.x),
                    std::sqrt(col.y),
                    std::sqrt(col.z));
            std::cout
                << int(255.99 * col.x) << " "
                << int(255.99 * col.y) << " "
                << int(255.99 * col.z) << " "
                << "\n";
        }
    }
        


    // Convert command-line arguments to vector of strings TODO
    // Declare variables for parsed command line TODO
    // Process command-line arguments TODO
    // Initialize pbrt TODO
    // Parse provided scene description files TODO
    // Render the scene TODO
    // Clean up after rendering the scene TODO
}
