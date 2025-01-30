#include "primitive.h"

namespace pbrt{

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
}
