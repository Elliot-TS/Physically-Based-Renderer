#pragma once
#include "interaction.h"
#include "samplers.h"
#include "ray.h"

namespace pbrt {
    class Material {
        public:
            Sampler *sampler;

            virtual bool scatter(
                    const Ray& r_in, 
                    const SurfaceInteraction& si, 
                    Vector3f& attenuation,
                    Ray& scattered) = 0;
    };


    class Lambertian : public Material {
        public:
            Vector3f albedo;
            Sampler *sampler;

            Lambertian(const Vector3f &albedo, Sampler *sampler) : albedo(albedo), sampler(sampler) {}
            virtual bool scatter(
                    const Ray& r_in, 
                    const SurfaceInteraction& si, 
                    Vector3f& attenuation,
                    Ray& scattered)
            {
                Point3f target = 
                    si.point + 
                    Vector3f(si.normal) +
                    sampler->s_Vector3f();

                scattered = Ray(si.point, target - si.point);
                attenuation = albedo;
                return true;
            }
    };
}
