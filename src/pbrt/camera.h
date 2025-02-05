#pragma once
#include "pbrt/ray.h"

namespace pbrt{

    /**
     * Define Camera Class
    **/
    class Camera{
        public:
            Point3f origin;
            Point3f lower_left_corner;
            Vector3f horizontal;
            Vector3f vertical;

            Camera() {
                lower_left_corner = Point3f(-2,-1,-1);
                horizontal = Vector3f(4,0,0);
                vertical = Vector3f(0,2,0);
                origin = Point3f(0,0,0);
            }
            Camera(Point3f origin, Vector3f direction, Float fov=67, Float aspectRatio = 2) : origin(origin) {
                // Find the up vector
                Vector3f up(0,1,0);
                if (direction.x == 0 && 
                    direction.y == 0)
                { 
                    if (direction.z == 0)
                        direction = Vector3f(0,0,-1);
                    else up = Vector3f(0,0,1);
                }
                
                // Find the horizontal and vertical vectors
                direction = Normalize(direction);
                horizontal = Cross(direction, up);
                vertical = Cross(horizontal, direction);

                // Find the field of view and aspectRatio
                Float height = std::atan(fov);
                horizontal *= height * aspectRatio;
                vertical *= height;

                // Find the lower left corner
                lower_left_corner = origin + direction - horizontal - vertical;

                // Scale the horizontal and vertical vectors
                horizontal *= 2;
                vertical *= 2;
            }

            inline Ray get_ray(const Float u, const Float v) const {
                return Ray(
                        origin,
                        lower_left_corner + u*horizontal + v*vertical - origin);
            }
    };

}
