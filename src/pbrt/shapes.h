#pragma once
#include <optional>
#include <vector>
#include <memory>
#include "pbrt/ray.h"
#include "pbrt/interaction.h"
#include "pbrt/util/buffercache.h"

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

        virtual Bounds3f Bounds() const = 0;
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

        Bounds3f Bounds() const {
            Vector3f diagonal(radius, radius, radius);
            return Bounds3f
                (center - diagonal,
                 center + diagonal);
        }
        
        // Members
        Float radius;
        Point3f center;
};


class TriangleMesh {
    public:
        int nTriangles, nVertices;
        const unsigned int *vertexIndeces = nullptr;
        const Point3f *vertexPositions = nullptr;

        TriangleMesh(
                std::vector<Point3f> positions,
                std::vector<unsigned int> indeces)
            : nTriangles(indeces.size()/3), nVertices(positions.size())
        {
            // TODO: Transform positions from object space to render space
            // Constructor should take a render space transform object
            vertexPositions = point3BufferCache->LookupOrAdd(positions);

            vertexIndeces = uintBufferCache->LookupOrAdd(indeces);
        }
};

struct TriangleIntersection {
    Float b0, b1, b2;
    Float t;
};

class Triangle : public Shape {
    public:
        Triangle (int meshIndex, int triIndex)
            : meshIndex(meshIndex), triIndex(triIndex) {}
        Bounds3f Bounds() const;

        static std::vector<std::unique_ptr<Shape>> CreateTriangles(const TriangleMesh *mesh);

        static SurfaceInteraction InteractionFromIntersection
            (const TriangleMesh *mesh, int triIndex,
             TriangleIntersection ti, Float time, Vector3f wo) 
        {
            const unsigned int *v = &mesh->vertexIndeces[3 * triIndex];
            Point3f p0 = mesh->vertexPositions[v[0]],
                    p1 = mesh->vertexPositions[v[1]],
                    p2 = mesh->vertexPositions[v[2]];

            // Find the hit point
            Point3f pHit = ti.b0*p0 + ti.b1*p1 + ti.b2*p2;

            // Find the normal TODO: Do better
            Normal3f n = Normal3f(Normalize(Cross((p1-p0), (p2-p0))));

            return SurfaceInteraction(pHit, n);
        }

        std::optional<ShapeIntersection> Intersect
            (const Ray &ray, Float tMax=Infinity) const;
        virtual bool IntersectP
            (const Ray &ray, Float tMax=Infinity) const;

        static void Init();

    private:
        int meshIndex = -1, triIndex = -1;
        static std::vector<const TriangleMesh *> *allMeshes;

        const TriangleMesh* GetMesh() const {
            return (*allMeshes)[meshIndex];
        }
        Float Area() const {
            const TriangleMesh *mesh = GetMesh();
            const unsigned int *v = &mesh->vertexIndeces[3 * triIndex];
            Point3f p0 = mesh->vertexPositions[v[0]],
                    p1 = mesh->vertexPositions[v[1]],
                    p2 = mesh->vertexPositions[v[2]];
            return 0.5f * Length(Cross(p1-p0, p2-p0));
        }

};
}

