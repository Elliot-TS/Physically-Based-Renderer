#include "aggregates.h"
#include <algorithm>
#include <atomic>

namespace pbrt {
// Simple Aggregate
std::optional<ShapeIntersection> SimpleAggregate::Intersect(
    const Ray &ray, Float tMax
) const
{
  Float closest = tMax;
  std::optional<ShapeIntersection> hit;
  for (int i = 0; i < numPrimitives; i++) {
    auto si = primitives[i]->Intersect(ray, closest);
    if (bool(si)) {
      hit = si;
      closest = si->tHit;
    }
  }
  return hit;
}
bool SimpleAggregate::IntersectP(const Ray &ray, Float tMax)
    const
{
  Float closest = tMax;
  for (int i = 0; i < numPrimitives; i++) {
    auto si = primitives[i]->IntersectP(ray, closest);
    if (bool(si)) return true;
  }
  return false;
}


// BVHAggregate
BVHAggregate::BVHAggregate(
    std::vector<Primitive *> prims, int maxPrimsInNode,
    SplitMethod splitMethod
)
    : maxPrimsInNode(std::min(255, maxPrimsInNode)),
      primitives(std::move(prims)),
      splitMethod(splitMethod)
{
  const int ps = 51;

  // Initialize bvhPrimitives array for primitives
  std::vector<BVHPrimitive> bvhPrimitives;
  for (size_t i = 0; i < primitives.size(); ++i) {
    bvhPrimitives.push_back(
        BVHPrimitive(i, primitives[i]->Bounds())
    );
  }

  // Build BVH for primitives using bvhPrimitives
  std::vector<Primitive *> orderedPrims(primitives.size());
  BVHBuildNode *root;

  std::atomic<int> totalNodes {0};
  if (splitMethod == SplitMethod::HLBVH) {
    // root = buildHLBVH(bvhPrimitives, &totalNodes,
    // orderedPrims);
    ERROR("SplitMethod::HLBVH is not yet implemented.");
  }
  else {
    std::atomic<int> orderedPrimsOffset {0};
    root = buildRecursive(
        bvhPrimitives, &totalNodes, &orderedPrimsOffset,
        orderedPrims
    );
  }

  primitives.swap(orderedPrims);

  // Convert BVH into compact representation in node array
  // bvhPrimitives.resize(0);  // free memory from bvhPrimitives
  // bvhPrimitives.shrink_to_fit();
  nodes = new LinearBVHNode[totalNodes];
  int offset = 0;
  flattenBVH(root, &offset);

  // TODO: Free BVHBuildNodes
}

inline void BVHAggregate::buildLeafNode(
    std::span<BVHPrimitive> &bvhPrimitives,
    std::atomic<int> *totalNodes,
    std::atomic<int> *orderedPrimsOffset,
    std::vector<Primitive *> &orderedPrims,
    Bounds3f bounds,
    BVHBuildNode *node
)
{
  int firstPrimOffset =
      orderedPrimsOffset->fetch_add(bvhPrimitives.size());
  for (size_t i = 0; i < bvhPrimitives.size(); ++i) {
    int index = bvhPrimitives[i].primitiveIndex;
    orderedPrims[firstPrimOffset + i] = primitives[index];
  }
  node->InitLeaf(firstPrimOffset, bvhPrimitives.size(), bounds);
}

BVHBuildNode *BVHAggregate::buildRecursive(
    std::span<BVHPrimitive> bvhPrimitives,
    std::atomic<int> *totalNodes,
    std::atomic<int> *orderedPrimsOffset,
    std::vector<Primitive *> &orderedPrims
)
{
  // Create a node to add
  BVHBuildNode *node = new BVHBuildNode;
  nodesToFree.push_back(node);
  ++*totalNodes;

  // Compute bounds for all primitives
  Bounds3f bounds;
  for (const BVHPrimitive &prim : bvhPrimitives) {
    bounds = Union(bounds, prim.bounds);
  }

  // Create a leaf node
  if (bounds.SurfaceArea() == 0 || bvhPrimitives.size() == 1) {
    buildLeafNode(
        bvhPrimitives, totalNodes, orderedPrimsOffset,
        orderedPrims, bounds, node
    );
    return node;
  }
  // Recursively create interior nodes
  else {
    // Compute bound of primitive centroids
    // and choose split dimension
    Bounds3f centroidBounds;
    for (const auto &prim : bvhPrimitives)
      centroidBounds = Union(centroidBounds, prim.Centroid());
    int dim = centroidBounds.MaxDimension();

    // If all the centroids occupy the same location,
    // create a leaf node
    if (centroidBounds.pMax[dim] == centroidBounds.pMin[dim]) {
      buildLeafNode(
          bvhPrimitives, totalNodes, orderedPrimsOffset,
          orderedPrims, bounds, node
      );
      return node;
    }
    // Partition primitives into two sets and build children
    else {
      int mid = bvhPrimitives.size() / 2;
      // Partition primitives based on splitMethod
      switch (splitMethod) {
        // Partition along middle of centroidBounds
        case SplitMethod::Middle: {
          if (partitionMiddle(
                  bvhPrimitives, centroidBounds, dim, &mid
              ))
            break;
        }
        // Partition primitives into equally sized subsets
        case SplitMethod::EqualCounts: {
          partitionEqualCounts(bvhPrimitives, dim, &mid);
          break;
        }
        case SplitMethod::SAH: {
          if (partitionSAH(
                  bvhPrimitives, centroidBounds, dim, &mid,
                  bounds, totalNodes, orderedPrimsOffset,
                  orderedPrims, bounds, node
              ))
          {
            return node;
          }
          break;
        }
        default:
          ERROR("The BVH Split Method is not implemented.");
      }

      // Recursively build BVHs for children
      // TODO: Make paralell if bvhPrimitives.size() > 128*1024
      BVHBuildNode *children[2];
      children[0] = buildRecursive(
          bvhPrimitives.subspan(0, mid), totalNodes,
          orderedPrimsOffset, orderedPrims
      );
      children[1] = buildRecursive(
          bvhPrimitives.subspan(mid), totalNodes,
          orderedPrimsOffset, orderedPrims
      );

      // Create the interior node
      node->InitInterior(dim, children[0], children[1]);
      return node;
    }
  }
}

inline bool BVHAggregate::partitionSAH(
    std::span<BVHPrimitive> &bvhPrimitives,
    Bounds3f centroidBounds, int dim, int *mid,
    Bounds3f primitiveBounds, std::atomic<int> *totalNodes,
    std::atomic<int> *orderedPrimsOffset,
    std::vector<Primitive *> &orderedPrims, Bounds3f bounds,
    BVHBuildNode *node
)
{
  // TODO: Currently working on 03/31 10:51
  if (bvhPrimitives.size() <= 2) {
    // Partition primitives into equally sized subsets
    partitionEqualCounts(bvhPrimitives, dim, mid);
  }
  else {
    // Allocate BVHSplitBucket for SAH partition buckets.
    constexpr int nBuckets = 12;
    BVHSplitBucket buckets[nBuckets];

    // Initialize BVHSplitBucket for SAH partition buckets.
    for (auto prim : bvhPrimitives) {
      // Determine which bucket index the primitive lies in
      // We multiply by 0.999 so bucket is the range
      // [0,nBuckets-1]
      int bucketIndex =
          (nBuckets * 0.999) *
          centroidBounds.Offset(prim.Centroid())[dim];
      buckets[bucketIndex].bounds =
          Union(buckets[bucketIndex].bounds, prim.bounds);
      buckets[bucketIndex].count++;
    }

    // Compute cost for split that minimizes SAH metric.
    // Cost of split A|B =
    //    1/2 prob_A * A.count + prob_B * B.count
    // TODO: See if you can do this in one pass
    int costs[nBuckets - 1] = {};
    int countBelow = 0;
    Bounds3f boundsBelow;
    for (int i = 0; i < nBuckets - 1; ++i) {
      countBelow += buckets[i].count;
      boundsBelow = Union(boundsBelow, buckets[i].bounds);
      costs[i] += boundsBelow.SurfaceArea() * countBelow;
    }

    int countAbove = 0;
    Bounds3f boundsAbove;
    for (int i = nBuckets - 2; i >= 0; --i) {
      countAbove += buckets[i].count;
      boundsAbove = Union(boundsAbove, buckets[i].bounds);
      costs[i] += boundsAbove.SurfaceArea() * countAbove;
    }

    // Find which split bucket would lead to the least cost
    int minCostBucket = -1;
    Float minCost = Infinity;
    for (int i = 0; i < nBuckets - 1; ++i) {
      if (costs[i] < minCost) {
        minCostBucket = i;
        minCost = costs[i];
      }
    }

    // Either create leaf or split primitives at selected SAH
    // bucket.
    Float leafCost = bvhPrimitives.size();
    minCost = 1 / 2.f + minCost / primitiveBounds.SurfaceArea();
    if (bvhPrimitives.size() > maxPrimsInNode ||
        minCost < leafCost)
    {
      auto midIter = std::partition(
          bvhPrimitives.begin(), bvhPrimitives.end(),
          [=](const BVHPrimitive &bp) {
            int bucketIndex =
                (nBuckets * 0.999) *
                centroidBounds.Offset(bp.Centroid())[dim];
            return bucketIndex <= minCostBucket;
          }
      );
      *mid = midIter - bvhPrimitives.begin();
    }
    else {
      buildLeafNode(
          bvhPrimitives, totalNodes, orderedPrimsOffset,
          orderedPrims, bounds, node
      );
      return true;
    }
  }
  return false;
}

inline bool BVHAggregate::partitionMiddle(
    std::span<BVHPrimitive> &bvhPrimitives,
    Bounds3f centroidBounds, int dim, int *mid
)
{
  Float pmid =
      (centroidBounds.pMin[dim] + centroidBounds.pMax[dim]) / 2;
  auto midIter = std::partition(
      bvhPrimitives.begin(), bvhPrimitives.end(),
      [dim, pmid](const BVHPrimitive &pi) {
        return pi.Centroid()[dim] < pmid;
      }
  );
  *mid = midIter - bvhPrimitives.begin();
  return midIter != bvhPrimitives.begin() &&
         midIter != bvhPrimitives.end();
}
inline void BVHAggregate::partitionEqualCounts(
    std::span<BVHPrimitive> &bvhPrimitives, int dim, int *mid
)
{
  *mid = bvhPrimitives.size() / 2;
  std::nth_element(
      bvhPrimitives.begin(), bvhPrimitives.begin() + *mid,
      bvhPrimitives.end(),
      [dim](const BVHPrimitive &a, const BVHPrimitive &b) {
        return a.Centroid()[dim] < b.Centroid()[dim];
      }
  );
}

int BVHAggregate::flattenBVH(BVHBuildNode *node, int *offset)
{
  // TODO: Would it be faster without recursion?
  int nodeOffset = *offset;
  LinearBVHNode *linearNode = &nodes[nodeOffset];
  linearNode->bounds = node->bounds;
  ++(*offset);
  if (node->nPrimitives > 0) {
    // Create leaf node
    linearNode->primitivesOffset = node->firstPrimOffset;
    linearNode->nPrimitives = node->nPrimitives;
  }
  else {
    // Recursively create interior flattened BVH node
    linearNode->axis = node->splitAxis;
    linearNode->nPrimitives = 0;
    flattenBVH(node->children[0], offset);
    linearNode->secondChildOffset =
        flattenBVH(node->children[1], offset);
  }
  return nodeOffset;
}

std::optional<ShapeIntersection> BVHAggregate::Intersect(
    const Ray &ray, Float tMax
) const
{
  std::optional<ShapeIntersection> si;
  Vector3f invDir(
      1 / ray.direction.x,
      1 / ray.direction.y,
      1 / ray.direction.z
  );
  int dirIsNeg[3] = {
      int(invDir.x < 0), int(invDir.y < 0), int(invDir.z < 0)
  };

  int numParents = 0;

  // Follow the ray through BVH nodes to find intersections
  int toVisitOffset = 0, currentNodeIndex = 0;
  int nodesToVisit[64];  // stack of nodes to visit
  while (true) {
    const LinearBVHNode *node = &nodes[currentNodeIndex];
    // If we intersect the current bounding box
    if (node->bounds.IntersectP(
            ray.origin, ray.direction, tMax, invDir, dirIsNeg
        ))
    {
      ++numParents;
      // Leaf Node
      if (node->nPrimitives > 0) {
        // Intersect ray with primitives in leaf node
        for (int i = 0; i < node->nPrimitives; ++i) {
          std::optional<ShapeIntersection> primSi =
              primitives[node->primitivesOffset + i]->Intersect(
                  ray, tMax
              );
          if (primSi) {
            si = primSi;
            tMax = si->tHit;
          }
        }

        // Continue visiting remaining nodes
        if (toVisitOffset == 0) break;
        currentNodeIndex = nodesToVisit[--toVisitOffset];
      }
      // Interior Node
      else {
        // Decide whether to visit first or second child first
        if (dirIsNeg[node->axis]) {  // second child first
          nodesToVisit[toVisitOffset++] = currentNodeIndex + 1;
          currentNodeIndex = node->secondChildOffset;
        }
        else {  // first child first
          nodesToVisit[toVisitOffset++] =
              node->secondChildOffset;
          ++currentNodeIndex;
        }
      }
    }
    else {
      if (toVisitOffset == 0) break;
      currentNodeIndex = nodesToVisit[--toVisitOffset];
    }
  }

  return si;
}
bool BVHAggregate::IntersectP(const Ray &ray, Float tMax) const
{
  Vector3f invDir(
      1 / ray.direction.x,
      1 / ray.direction.y,
      1 / ray.direction.z
  );
  int dirIsNeg[3] = {
      int(invDir.x < 0), int(invDir.y < 0), int(invDir.z < 0)
  };

  // Follow the ray through BVH nodes to find intersections
  int toVisitOffset = 0, currentNodeIndex = 0;
  int nodesToVisit[64];  // stack of nodes to visit
  while (true) {
    const LinearBVHNode *node = &nodes[currentNodeIndex];
    // If we intersect the current bounding box
    if (node->bounds.IntersectP(
            ray.origin, ray.direction, tMax, invDir, dirIsNeg
        ))
    {
      // Leaf Node
      if (node->nPrimitives > 0) {
        // Intersect ray with primitives in leaf node
        for (int i = 0; i < node->nPrimitives; ++i) {
          bool hit = primitives[node->primitivesOffset + i]
                         ->IntersectP(ray, tMax);
          if (hit) return true;
        }

        // Continue visiting remaining nodes
        if (toVisitOffset == 0) break;
        currentNodeIndex = nodesToVisit[--toVisitOffset];
      }
      // Interior Node
      else {
        // Decide whether to visit first or second child first
        if (dirIsNeg[node->axis]) {  // second child first
          nodesToVisit[toVisitOffset++] = currentNodeIndex + 1;
          currentNodeIndex = node->secondChildOffset;
        }
        else {  // first child first
          nodesToVisit[toVisitOffset++] =
              node->secondChildOffset;
          ++currentNodeIndex;
        }
      }
    }
    else {
      if (toVisitOffset == 0) break;
      currentNodeIndex = nodesToVisit[--toVisitOffset];
    }
  }

  return false;
}
}  // namespace pbrt
