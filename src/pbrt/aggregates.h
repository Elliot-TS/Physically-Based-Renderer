#pragma once
#include "primitive.h"

namespace pbrt {

struct BVHPrimitive {
  size_t primitiveIndex;
  Bounds3f bounds;

  BVHPrimitive(size_t primitiveIndex, const Bounds3f &bounds)
      : primitiveIndex(primitiveIndex), bounds(bounds)
  {}
  BVHPrimitive() = default;

  Point3f Centroid() const
  {
    return .5f * bounds.pMin + 0.5f * bounds.pMax;
  }
};
struct BVHBuildNode {
  void InitLeaf(int first, int n, const Bounds3f &b)
  {
    firstPrimOffset = first;
    nPrimitives = n;
    bounds = b;
    children[0] = children[1] = nullptr;
  }
  void InitInterior(
      int axis, BVHBuildNode *c0, BVHBuildNode *c1
  )
  {
    children[0] = c0;
    children[1] = c1;
    bounds = Union(c0->bounds, c1->bounds);
    splitAxis = axis;
    nPrimitives = 0;
  }
  Bounds3f bounds;
  BVHBuildNode *children[2];
  int splitAxis, firstPrimOffset, nPrimitives;
};
struct alignas(32) LinearBVHNode {
  Bounds3f bounds;
  union {
    int primitivesOffset;   // leaf
    int secondChildOffset;  // interior
  };
  uint16_t nPrimitives;
  uint8_t axis;
};
struct BVHSplitBucket {
  int count = 0;
  Bounds3f bounds;
};

class SimpleAggregate : public Primitive {
 public:
  // TODO: Change to vector, for consistency with BVHAggregate
  Primitive **primitives;  // List of pointers to primitives
  int numPrimitives;

  SimpleAggregate(Primitive **primitives, int numPrimitives)
      : primitives(primitives), numPrimitives(numPrimitives)
  {}

  std::optional<ShapeIntersection> Intersect(
      const Ray &ray, Float tMax = Infinity
  ) const;
  bool IntersectP(const Ray &ray, Float tMax = Infinity) const;

  Bounds3f Bounds() const
  {
    Bounds3f totalBounds;
    for (int i = 0; i < numPrimitives; ++i) {
      totalBounds = Union(totalBounds, primitives[i]->Bounds());
    }
    return totalBounds;
  }
};

class BVHAggregate : public Primitive {
 private:
  std::vector<BVHBuildNode *> nodesToFree;

  void buildLeafNode(
      std::span<BVHPrimitive> &bvhPrimitives,
      std::atomic<int> *totalNodes,
      std::atomic<int> *orderedPrimsOffset,
      std::vector<Primitive *> &orderedPrims,
      Bounds3f bounds,
      BVHBuildNode *node
  );

  bool partitionMiddle(
      std::span<BVHPrimitive> &bvhPrimitives,
      Bounds3f centroidBounds, int dim, int *mid
  );
  // Returns whether a leaf node was created
  bool partitionSAH(
      std::span<BVHPrimitive> &bvhPrimitives,
      Bounds3f centroidBounds, int dim, int *mid,
      Bounds3f primitiveBounds, std::atomic<int> *totalNodes,
      std::atomic<int> *orderedPrimsOffset,
      std::vector<Primitive *> &orderedPrims, Bounds3f bounds,
      BVHBuildNode *node
  );
  void partitionEqualCounts(
      std::span<BVHPrimitive> &bvhPrimitives, int dim, int *mid
  );

 public:
  // Algorithms for partitioning primitives
  enum class SplitMethod { SAH, HLBVH, Middle, EqualCounts };

  // Public members
  int maxPrimsInNode;
  std::vector<Primitive *> primitives;
  SplitMethod splitMethod;
  LinearBVHNode *nodes = nullptr;  // probably should be private

  BVHAggregate(
      std::vector<Primitive *> prims, int maxPrimsInNode,
      SplitMethod splitMethod
  );
  ~BVHAggregate()
  {
    delete[] nodes;
    for (auto node : nodesToFree) {
      delete node;
    }
  }

  BVHBuildNode *buildRecursive(
      std::span<BVHPrimitive> bvhPrimitives,
      std::atomic<int> *totalNodes,
      std::atomic<int> *orderedPrimsOffset,
      std::vector<Primitive *> &orderedPrims
  );

  int flattenBVH(BVHBuildNode *node, int *offset);

  Bounds3f Bounds() const { return nodes[0].bounds; }

  std::optional<ShapeIntersection> Intersect(
      const Ray &ray, Float tMax
  ) const;
  std::optional<ShapeIntersection> Intersect(
      const Ray &ray, Float tMax, double *trackBVHLayers
  ) const;
  bool IntersectP(const Ray &ray, Float tMax) const;
};

}  // namespace pbrt
