#pragma once
#include <optional>
#include "pbrt/math/vecmath.h"
#include "pbrt/primitive.h"
#include "pbrt/shapes.h"
#include "pbrt/light.h"
#include "pbrt/samplers.h"
#include "pbrt/camera.h"

namespace pbrt {

    class Integrator {
        public:
            // Methods
            Integrator(Primitive *aggregate, std::vector<Light> lights)
                : aggregate(aggregate), lights(lights)
            {
                // TODO: Integrator Constructor Implementation
            }

            virtual void Render() = 0;
            std::optional<ShapeIntersection> Intersect(const Ray &ray, Float tMax) const;
            bool IntersectP(const Ray &ray, Float tMax) const;

            // Members TODO
            Primitive *aggregate;
            std::vector<Light> lights;
            std::vector<Light> infiniteLights;
        protected:
            // Methods TODO
    };

    class ImageTileIntegrator : public Integrator {
        public:
            Camera *camera;
            Sampler *samplerPrototype;

            ImageTileIntegrator(
                    Camera *camera,
                    Sampler *sampler,
                    Primitive *aggregate,
                    std::vector<Light> lights):
                Integrator(aggregate, lights), camera(camera), samplerPrototype(sampler) {}

            void Render();
            //virtual void EvaluatePixelSample(Point2f pPixel, int sampleIndex, Sampler *sampler) = 0;
    };
}
