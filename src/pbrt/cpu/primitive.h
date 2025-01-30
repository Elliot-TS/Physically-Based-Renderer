// Based off PBRT and RTW, but also my own code
#pragma once
#include "../ray.h"
#include "../shapes.h"

namespace pbrt {

    class Primitive {
        public:
            // Bounds3f Bounds() const;
            virtual std::optional<ShapeIntersection> Intersect
                (const Ray &ray, Float tMax = Infinity) const = 0;
            virtual bool IntersectP
                (const Ray &ray, Float tMax = Infinity) const = 0;

            // Public Members
            // Material material;
            // Light areaLight;
    };


    class GeometricPrimitive : public Primitive {
        public:
            Shape *shape;

            GeometricPrimitive(Shape *shape): shape(shape) {}

            std::optional<ShapeIntersection> Intersect
                (const Ray &ray, Float tMax) const {
                    return shape->Intersect(ray, tMax);
                }
            bool IntersectP (const Ray &ray, Float tMax) const {
                return shape->IntersectP(ray, tMax);
            }
    };

}
