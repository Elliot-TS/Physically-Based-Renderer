// Based off PBRT and RTW, but also my own code
#pragma once
#include "pbrt/ray.h"
#include "pbrt/shapes.h"
#include "pbrt/materials.h"

namespace pbrt {

    class Primitive {
        public:
            // Bounds3f Bounds() const;
            virtual std::optional<ShapeIntersection> Intersect
                (const Ray &ray, Float tMax = Infinity) const = 0;
            virtual bool IntersectP
                (const Ray &ray, Float tMax = Infinity) const = 0;

            // Public Members
            Material *material;
            // Light areaLight;
    };


    class GeometricPrimitive : public Primitive {
        public:
            Shape *shape;
            Material *material;

            GeometricPrimitive(Shape *shape, Material *material): 
                shape(shape), material(material) {}

            std::optional<ShapeIntersection> Intersect
                (const Ray &ray, Float tMax) const {
                    std::optional<ShapeIntersection> si = shape->Intersect(ray, tMax);
                    if (!si) return {};
                    si->interaction.material = material;
                    return si;
                }
            bool IntersectP (const Ray &ray, Float tMax) const {
                return shape->IntersectP(ray, tMax);
            }
    };

    class SimpleAggregate : public Primitive {
        public:
            Primitive **primitives; // List of pointers to primitives
            int numPrimitives;

            SimpleAggregate
                (Primitive **primitives, int numPrimitives):
                    primitives(primitives),
                    numPrimitives(numPrimitives) {}

            std::optional<ShapeIntersection> Intersect
                (const Ray &ray, Float tMax = Infinity) const;
            bool IntersectP
                (const Ray &ray, Float tMax = Infinity) const;
    };
}
