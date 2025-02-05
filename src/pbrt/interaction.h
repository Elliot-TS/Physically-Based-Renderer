#pragma once
#include "pbrt.h"
#include "pbrt/math/vecmath.h"

namespace pbrt {
class SurfaceInteraction {
public:
    // Public Methods
    SurfaceInteraction() {};
    SurfaceInteraction(Point3f point, Normal3f normal):
        point(point), normal(normal) {};

    // Public Members
    Point3f point;
    Normal3f normal;
    Material *material;
};
}
