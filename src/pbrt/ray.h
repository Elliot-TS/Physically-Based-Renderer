// Mostly copied from PBRT
#pragma once
#include "pbrt/math/vecmath.h"

namespace pbrt {
// Ray Definition
// TODO: Implement Medium
class Ray {
 public:
  // Ray Public Methods
  bool HasNaN() const
  {
    return (origin.HasNaN() || direction.HasNaN());
  }

  friend std::ostream& operator<<(
      std::ostream& os, const Ray ray
  )
  {
    os << "["
       << " o: " << ray.origin << " d: " << ray.direction
       << " time: "
       << ray.time
       // << " medium: " << ray.medium
       << " ]";
    return os;
  }

  Point3f operator()(Float t) const
  {
    return origin + direction * t;
  }

  Ray() = default;
  Ray(Point3f origin,
      Vector3f direction,
      Float time = 0.f /*, Medium medium = nullptr*/)
      : origin(origin),
        direction(direction),
        time(time) /*, medium(medium)*/
  {}

  // Ray Public Members
  Point3f origin;
  Vector3f direction;
  Float time = 0;
  // Medium medium = nullptr;
};

// RayDifferential Definition
class RayDifferential : public Ray {
 public:
  // RayDifferential Public Methods
  RayDifferential() = default;

  RayDifferential(
      Point3f o, Vector3f d, Float time = 0.f
      // Medium medium = nullptr
  )
      : Ray(o, d, time /*, medium*/)
  {}


  explicit RayDifferential(const Ray& ray): Ray(ray) {}

  void ScaleDifferentials(Float s)
  {
    rxOrigin = origin + (rxOrigin - origin) * s;
    ryOrigin = origin + (ryOrigin - origin) * s;
    rxDirection = direction + (rxDirection - direction) * s;
    ryDirection = direction + (ryDirection - direction) * s;
  }


  bool HasNaN() const
  {
    return Ray::HasNaN() ||
           (hasDifferentials &&
            (rxOrigin.HasNaN() || ryOrigin.HasNaN() ||
             rxDirection.HasNaN() || ryDirection.HasNaN()));
  }

  // RayDifferential Public Members
  bool hasDifferentials = false;
  Point3f rxOrigin, ryOrigin;
  Vector3f rxDirection, ryDirection;
};
}  // namespace pbrt
