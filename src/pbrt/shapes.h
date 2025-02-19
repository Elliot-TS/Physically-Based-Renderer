#pragma once
#include <optional>
#include <vector>
#include "pbrt/ray.h"
#include "pbrt/interaction.h"

namespace pbrt {
    
    struct ShapeIntersection {
        SurfaceInteraction interaction;
        Float tHit;

        ShapeIntersection(SurfaceInteraction interaction, Float tHit):
            interaction(interaction), tHit(tHit) {}
        ShapeIntersection() {}
    };

    class Shape {
        public:
            virtual std::optional<ShapeIntersection> Intersect
                (const Ray &ray, Float tMax=Infinity) const = 0;
            virtual bool IntersectP
                (const Ray &ray, Float tMax=Infinity) const = 0;
    };


    class Sphere : public Shape { 
        private:
            Float tMin = 0.001;
        public:
            // Methods
            Sphere() {}
            Sphere(Point3f center, Float radius) : 
                center(center), radius(radius) {};

            std::optional<ShapeIntersection> Intersect
                (const Ray &ray, Float tMax=Infinity) const;
            inline bool IntersectP
                (const Ray &ray, Float tMax=Infinity) const;
            
            // Members
            Float radius;
            Point3f center;
    };


    class TriangleMesh : Shape {
        public:
            int nTriangles, nVertices;
            const int *vertexIndeces = nullptr;
            const Vector3f *vertexPositions = nullptr;

            TriangleMesh(
                    std::vector<Point3f> positions,
                    std::vector<int> indeces)
                : nTriangles(indeces.size()/3), nVertices(positions.size())
            {
            }
    };
}
