#include "aggregates.h"
#include <atomic>
#include <algorithm>

namespace pbrt {
    // Simple Aggregate
    std::optional<ShapeIntersection> SimpleAggregate::Intersect
        (const Ray &ray, Float tMax) const
        {
            Float closest = tMax;
            std::optional<ShapeIntersection> hit;
            for (int i = 0; i < numPrimitives; i++) {
                auto si = primitives[i]->Intersect
                    (ray, closest);
                if (bool(si)) {
                    hit = si;
                    closest = si->tHit;
                }
            }    
            return hit;
        }
    bool SimpleAggregate::IntersectP
        (const Ray &ray, Float tMax) const
        {
            Float closest = tMax;
            for (int i = 0; i < numPrimitives; i++) {
                auto si = primitives[i]->IntersectP
                    (ray, closest);
                if (bool(si))
                    return true;
            }
            return false;
        }


    // BVHAggregate
    BVHAggregate::BVHAggregate
        ( std::vector<Primitive*> primitives,
          int maxPrimsInNode,
          SplitMethod splitMethod
        ):   
            maxPrimsInNode(std::min(255, maxPrimsInNode)),
            primitives(std::move(primitives)),
            splitMethod(splitMethod)
    {
        // Initialize bvhPrimitives array for primitives
        std::vector<BVHPrimitive> bvhPrimitives;
        for (size_t i = 0; i < primitives.size(); ++i) {
            bvhPrimitives.push_back(
                    BVHPrimitive(i, primitives[i]->Bounds()));
        }

        // Build BVH for primitives using bvhPrimitices
        // Maybe TODO: Declare Allocators used for BVH construction
        std::vector<Primitive*> orderedPrims;
        BVHBuildNode *root;

        // Build BVH according to selected splitMethod
        std::atomic<int> totalNodes{0};
        if (splitMethod == SplitMethod::HLBVH) {
            ERROR("SplitMethod::HLBVH is not implemented.  Choose a different method for generating the BVH");
            //root = buildHLBVH(bhvPrimitives, &totalNodes, orderedPrims);
        }
        else {
            std::atomic<int> orderedPrimsOffset{0};
            root = buildRecursive(
                    std::span<BVHPrimitive>(bvhPrimitives),
                    &totalNodes, &orderedPrimsOffset, orderedPrims);
        }
        primitives.swap(orderedPrims);

        // Convert BHV into compact representation in nodes array
    }


    BVHBuildNode *BVHAggregate::buildRecursive(
            std::span<BVHPrimitive> bvhPrimitives,
            std::atomic<int> *totalNodes,
            std::atomic<int> *orderedPrimsOffset,
            std::vector<Primitive*> &orderedPrims) 
    {
        // TODO: Maybe use allocators
        BVHBuildNode *node = new BVHBuildNode;

        // Initialize BVHBuildNode for all primitive range
        ++(*totalNodes);
        Bounds3f bounds;

        auto createLeafNode = [&](BVHBuildNode *node){
            int firstPrimOffset = orderedPrimsOffset->fetch_add(bvhPrimitives.size());
            for (size_t i = 0; i < bvhPrimitives.size(); ++i) {
                int index = bvhPrimitives[i].primitiveIndex;
                orderedPrims[firstPrimOffset + i] = primitives[index];
            }
            node->InitLeaf(firstPrimOffset, bvhPrimitives.size(), bounds);
        };

        // Compute the bounds of all the primitives in the node
        for (const auto &prim: bvhPrimitives)
            bounds = Union(bounds, prim.bounds);

        // If there is only one primitive, create a leaf
        if (bounds.SurfaceArea() == 0 || bvhPrimitives.size() == 1) {
            // Create leaf node
            createLeafNode(node);
            return node;
        }
        else {
            // Compute bounds of primitive centroids and choose a split dimension
            Bounds3f centroidBounds;
            for (const auto &prim: bvhPrimitives)
                centroidBounds = Union(centroidBounds, prim.Centroid());
            int dim = centroidBounds.MaxDimension();

            if (centroidBounds.pMax[dim] == centroidBounds.pMin[dim]) {
                createLeafNode(node);
                return node;
            }
            else {
                int mid = bvhPrimitives.size() / 2;

                // TODO: Partition primitives based on splitMethod
                switch (splitMethod) {
                    case SplitMethod::Middle: {
                        Float pmid = (centroidBounds.pMin[dim] + centroidBounds.pMax[dim]) / 2;
                        auto midIter =
                            std::partition(
                                    bvhPrimitives.begin(),
                                    bvhPrimitives.end(),
                                    [dim, pmid](const BVHPrimitive &pi) {
                                        return pi.Centroid()[dim] < pmid;
                                    });
                        mid = midIter - bvhPrimitives.begin();
                        if (midIter != bvhPrimitives.begin() &&
                            midIter != bvhPrimitives.end()) 
                        { break; }
                    }
                    default:
                        ERROR("The BVH Split Method is not implemented.");
                }
                
                
                BVHBuildNode *children[2];
                // Recursively build child BVHs
                // TODO: Make this parallel if bvhPrimitives.size() > 128*1024
                children[0] = buildRecursive(bvhPrimitives.subspan(0, mid), totalNodes, orderedPrimsOffset, orderedPrims);
                children[1] = buildRecursive(bvhPrimitives.subspan(mid), totalNodes, orderedPrimsOffset, orderedPrims);

                node->InitInterior(dim, children[0], children[1]);
            }

            // Partition primitives into two sets and build children
        }

        return node;
    }
}
