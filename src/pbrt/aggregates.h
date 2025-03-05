#pragma once
#include "primitive.h"

namespace pbrt {


    struct BVHPrimitive {
        size_t primitiveIndex;
        Bounds3f bounds;

        BVHPrimitive(size_t primitiveIndex, const Bounds3f &bounds)
            : primitiveIndex(primitiveIndex), bounds(bounds) {}
        BVHPrimitive() = default;

        Point3f Centroid() const { return .5f * bounds.pMin + 0.5f * bounds.pMax; }
    };
    struct BVHBuildNode {
        void InitLeaf(int first, int n, const Bounds3f &b) {
            firstPrimOffset = first;
            nPrimitives = n;
            bounds = b;
            children[0] = children[1] = nullptr;
        }
        void InitInterior(int axis, BVHBuildNode *c0, BVHBuildNode *c1) {
            children[0] = c0;
            children[1] = c1;
            bounds = Union(c0->bounds, c1->bounds);
            splitAxis = splitAxis;
            nPrimitives = 0;
        }
        Bounds3f bounds;
        BVHBuildNode *children[2];
        int splitAxis, firstPrimOffset, nPrimitives;
    };
    struct alignas(32) LinearBVHNode {
        Bounds3f bounds;
        union {
            int primitivesOffset; // leaf
            int secondChildOffset; // interior
        };
        uint16_t nPrimitives;
        uint8_t axis;
    };


    class SimpleAggregate : public Primitive {
        public:
            // TODO: Change to vector, for consistency with BVHAggregate
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

            Bounds3f Bounds() const {
                Bounds3f totalBounds;
                for (int i = 0; i < numPrimitives; ++i) {
                    totalBounds = Union(totalBounds, primitives[i]->Bounds());
                }
                return totalBounds;
            }
    };


    class BVHAggregate: public Primitive {
        public:
            // Algorithms for partitioning primitives
            enum class SplitMethod { SAH, HLBVH, Middle, EqualCounts };

            // Public members
            int maxPrimsInNode;
            std::vector<Primitive*> primitives;
            SplitMethod splitMethod;
            LinearBVHNode *nodes = nullptr;

            BVHAggregate
                ( std::vector<Primitive*> primitives,
                  int maxPrimsInNode,
                  SplitMethod splitMethod);

            BVHBuildNode *buildRecursive(
                    std::span<BVHPrimitive> bvhPrimitives,
                    std::atomic<int> *totalNodes,
                    std::atomic<int> *orderedPrimsOffset,
                    std::vector<Primitive*> &orderedPrims);

            int flattenBVH(BVHBuildNode *node, int *offset);
            
            Bounds3f Bounds() const {
                return nodes[0].bounds;
            }

            std::optional<ShapeIntersection> Intersect
                (const Ray &ray, Float tMax) const;
            bool IntersectP
                (const Ray &ray, Float tMax) const;
    };

}
