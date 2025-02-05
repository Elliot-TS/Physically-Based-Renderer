// Adapted from RTW and PBRT
#include "pbrt/shapes.h"
namespace pbrt {
    // Sphere
    std::optional<ShapeIntersection> Sphere::Intersect
        (const Ray &ray, Float tMax) const
        {
            Vector3f oc = ray.origin - center;

            Float a = Dot(ray.direction, ray.direction);
            Float b = Dot(oc, ray.direction);
            Float c = Dot(oc, oc) - radius*radius;
            Float discriminant = b*b - a*c;

            if (discriminant > 0) {
                Float temp = (-b - sqrt(b*b - a*c)) / a;
                if (temp < tMax && temp > tMin) {
                    // Surface Interaction contains
                    //  point and normal
                    SurfaceInteraction surfIntrc;
                    surfIntrc.point = ray(temp);
                    surfIntrc.normal = Normal3f
                        ((surfIntrc.point - center) / radius);

                    // Shape Intersection additionally contains
                    //  tHit
                    ShapeIntersection shapeIntrs;
                    shapeIntrs.tHit = temp;
                    shapeIntrs.interaction = surfIntrc;

                    return shapeIntrs;
                }

                temp = (-b + sqrt(b*b - a*c)) / a;
                if (temp < tMax && temp > tMin) {
                    SurfaceInteraction surfIntrc;
                    surfIntrc.point = ray(temp);
                    surfIntrc.normal = Normal3f
                        ((surfIntrc.point - center) / radius);

                    ShapeIntersection shapeIntrs;
                    shapeIntrs.tHit = temp;
                    shapeIntrs.interaction = surfIntrc;

                    return shapeIntrs;
                }
            }

            return {};
        }

    inline bool Sphere::IntersectP
        (const Ray &ray, Float tMax) const {
            return bool(Intersect(ray, tMax));
        }
}
