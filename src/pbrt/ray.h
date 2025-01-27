#include "util/vecmath.h"

namespace pbrt{
// Ray Definition
// TODO: Implement Medium
class Ray {
  public:
    // Ray Public Methods
    bool HasNaN() const { return (origin.HasNaN() || direction.HasNaN()); }

    friend std::ostream& operator<<(std::ostream& os, const Ray ray) {
        os << "["
            << " o: " << ray.origin
            << " d: " << ray.direction
            << " time: " << ray.time
            // << " medium: " << ray.medium
            << " ]";
        return os;
    }
    
    Point3f operator()(Float t) const { return origin + direction * t; }

    Ray() = default;
    Ray(Point3f origin, Vector3f direction, Float time = 0.f/*, Medium medium = nullptr*/)
        : origin(origin), direction(direction), time(time)/*, medium(medium)*/ {}

    // Ray Public Members
    Point3f origin;
    Vector3f direction;
    Float time = 0;
    //Medium medium = nullptr;
};

}
