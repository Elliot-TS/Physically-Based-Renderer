// Based off PBRT and RTW, but also my own code
#pragma once
#include "pbrt/materials.h"
#include "pbrt/ray.h"
#include "pbrt/shapes.h"

namespace pbrt {

class Primitive {
 public:
  // Bounds3f Bounds() const;
  virtual std::optional<ShapeIntersection> Intersect(
      const Ray &ray, Float tMax = Infinity
  ) const = 0;

  // This overload is for debugging BVH Aggregates by visualizing
  // BVH bounding boxes by multiplying trackBVHLayers by .98 every
  // time it intersects a new bounding box virtual
  std::optional<ShapeIntersection> Intersect(
      const Ray &ray, Float tMax, double *trackBVHLayers
  ) const
  {
    return {};
  }

  virtual bool IntersectP(const Ray &ray, Float tMax = Infinity)
      const = 0;

  virtual Bounds3f Bounds() const = 0;

  // Public Members
  Material *material;
  // Light areaLight;
};


class GeometricPrimitive : public Primitive {
 public:
  Shape *shape;
  Material *material;

  GeometricPrimitive(Shape *shape, Material *material)
      : shape(shape), material(material)
  {}

  std::optional<ShapeIntersection> Intersect(
      const Ray &ray, Float tMax
  ) const
  {
    std::optional<ShapeIntersection> si =
        shape->Intersect(ray, tMax);
    if (!si) return {};
    si->interaction.material = material;
    return si;
  }
  bool IntersectP(const Ray &ray, Float tMax) const
  {
    return shape->IntersectP(ray, tMax);
  }
  Bounds3f Bounds() const { return shape->Bounds(); }
};
}  // namespace pbrt
