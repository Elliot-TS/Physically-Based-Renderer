#include <iostream>
#include "pbrt/ray.h"

using namespace pbrt;

bool hit_sphere(const Point3f center, Float radius, const Ray& r) {
    Vector3f oc = r.origin - center;
    Float a = Dot(r.direction, r.direction);
    Float b = 2.0 * Dot(oc, r.direction);
    Float c = Dot(oc, oc) - Sqr(radius);
    Float discriminant = b*b - 4*a*c;
    return (discriminant > 0);
}

Vector3f color(const Ray& r) {
    if (hit_sphere(Point3f(0,0,-1), 0.5, r))
        return Vector3f(1,0,0);

    Vector3f unit_direction = Normalize(r.direction);
    Float t = 0.5*(unit_direction.y + 1.0);
    return (1.0-t)  * Vector3f(1,1,1) 
            + t     * Vector3f(0.5, 0.7, 1.0);
}

int main(int arcg, char *argv[]) {
    // Ray tracing in a weekend
    int nx = 200;
    int ny = 100;
    std::cout << "P3\n" << nx << " " << ny << "\n255\n";

    Vector3f lower_left_corner(-2, -1, -1);
    Vector3f horizontal(4, 0, 0);
    Vector3f vertical(0, 2, 0);
    Point3f origin(0,0,0);

    for (int j = ny-1; j >= 0; j--) {
        for (int i = 0; i < nx; i++) {
            Float u = Float(i) / Float(nx);
            Float v = Float(j) / Float(ny);
            
            Ray r(
                origin, 
                lower_left_corner + u*horizontal + v*vertical
            );
            Vector3f col = color(r);

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
