#include "util/vecmath.h"

namespace pbrt{
// Ray Definition
// TODO: Implement Medium
class Ray {
  public:
    // Ray Public Methods
    bool HasNaN() const { return (o.HasNaN() || d.HasNaN()); }

    friend std::ostream& operator<<(std::ostream& os, const Ray ray) {
        os << "["
            << " o: " << ray.o
            << " d: " << ray.d
            << " time: " << ray.time
            // << " medium: " << ray.medium
            << " ]";
        return os;
    }
    
    Point3f operator()(Float t) const { return o + d * t; }

    Ray() = default;
    Ray(Point3f o, Vector3f d, Float time = 0.f/*, Medium medium = nullptr*/)
        : o(o), d(d), time(time)/*, medium(medium)*/ {}

    // Ray Public Members
    Point3f o;
    Vector3f d;
    Float time = 0;
    //Medium medium = nullptr;
};

}
