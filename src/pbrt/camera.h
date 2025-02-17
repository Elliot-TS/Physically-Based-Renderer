#pragma once
#include "pbrt/ray.h"
#include "pbrt/util/display.h"

namespace pbrt{

    class Film {
        public:
            const int width;
            const int height;
            Vector3f *imageSamples; // TODO: Consider Float
            DisplayWindow *display;

            // TODO: Something for saving the image to a file
            Film(const int width, const int height) 
                : width(width), height(height)
            { 
                imageSamples = new Vector3f[width*height];
                display = new DisplayWindow(width, height);
            }

            // Adds value to the average of the last numSamples samples
            void AddSample(Vector3f color, int x, int y, int numSamples) {
                double invSamp = 1.0 / (numSamples + 1.0);
                double sampRatio = double(numSamples) * invSamp;
                int index = y*width + x;
                Vector3f col = imageSamples[index]*sampRatio + color*invSamp;
                imageSamples[index] = col; 
                display->SetPixel(index, Color(col));
            }

            void UpdateDisplay() {
                display->UpdateImage(imageSamples, width*height);
            }
    };

    /**
     * Define Camera Class
    **/
    class Camera{
        public:
            Point3f origin;
            Point3f lower_left_corner;
            Vector3f horizontal;
            Vector3f vertical;
            Film *film; // I'm not understanding why I have to use a pointer

            // TODO: Find a more elegant way of passing display
            Camera(Film *film): film(film) {
                lower_left_corner = Point3f(-2,-1,-1);
                horizontal = Vector3f(4,0,0);
                vertical = Vector3f(0,2,0);
                origin = Point3f(0,0,0);
            }

            Camera(Film *film, Point3f origin, Vector3f direction, Float fov=67, Float aspectRatio = 2) : origin(origin), film(film) {
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

            Ray get_ray(const Float u, const Float v) const {
                return Ray(
                        origin,
                        lower_left_corner + u*horizontal + v*vertical - origin);
            }
    };

}
