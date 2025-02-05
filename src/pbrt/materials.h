#pragma once
#include "pbrt/interaction.h"
#include "pbrt/samplers.h"
#include "pbrt/ray.h"

namespace pbrt {
    class Material {
        public:
            Sampler *sampler;

            virtual bool scatter(
                    const Ray& r_in, 
                    const SurfaceInteraction& si, 
                    Vector3f& attenuation,
                    Ray& scattered) const = 0;
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
                    Ray& scattered) const
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


    class Metal : public Material {
        public:
            Vector3f albedo;
            Sampler *sampler;
            Float roughness;

            Metal(const Vector3f albedo, Float roughness, Sampler *sampler): 
                albedo(albedo), roughness(roughness), sampler(sampler) { roughness = roughness > 1 ? 1 : roughness; }
            virtual bool scatter(
                    const Ray& r_in, 
                    const SurfaceInteraction& si, 
                    Vector3f& attenuation,
                    Ray& scattered) const
            {
                Vector3f reflected = Reflect(Normalize(r_in.direction), si.normal);
                scattered = Ray(si.point, reflected + Vector3f(roughness*sampler->s_uSphere()));
                attenuation = albedo;
                return (Dot(scattered.direction, si.normal) > 0);
            }
    };
}
