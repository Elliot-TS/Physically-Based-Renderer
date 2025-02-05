#include "pbrt/integrator.h"

namespace pbrt {
Vector3f color(const Ray& r, const Primitive& shape, int bounces) {
    auto si = shape.Intersect(r);
    if (si) {
        Ray scattered;
        Vector3f attenuation;
        if (bounces > 0 && 
                si->interaction.material->scatter(
                    r, si->interaction, attenuation, scattered)
           )
        {
            return HorizontalProduct(
                    attenuation,
                    color(scattered, shape, bounces+1));
        }
        else return Vector3f(0,0,0);
    }

    Vector3f unit_direction= Normalize(r.direction);
    Float t = 0.5*(unit_direction.y + 1.0);
    return (1.0-t)  * Vector3f(1,1,1) 
            + t     * Vector3f(0.5, 0.7, 1.0);
}

    /** Integrator **/
    std::optional<ShapeIntersection> Integrator::Intersect(const Ray &ray,  Float tMax) const {
        if (aggregate != nullptr) {
            return aggregate->Intersect(ray, tMax);
        }
        else { return {}; }
    }
    bool Integrator::IntersectP(const Ray &ray, Float tMax) const {
        if (aggregate != nullptr) {
            return aggregate->IntersectP(ray, tMax);
        }
        else { return false; }
    }

    /** ImageTileIntegrator **/
    void ImageTileIntegrator::Render() {
        // TODO: Declare common variables for rendering image in tiles
        // TODO: Multi-threading
        int nx = 200;
        int ny = 100;
        Sampler *sampler = samplerPrototype; // We're going to use a ScratchBuffer later, that's why we're making this alias
        
        // Render image in waves
        int samplesPerPixel = 256;
        int waveStart = 0, waveEnd = samplesPerPixel, nextWaveSize = 1;
        std::cout << "P3\n" << nx << " " << ny << "\n255\n";
        while (waveStart < samplesPerPixel)
        {
            // Render current wave's image tiles in parallel
            // TODO: Use mutli-threading via ParallelFor2D() function.  For now, assume there is just one tile.

            for (int j = ny-1; j >= 0; j--) {
                for (int i = 0; i < nx; i++) {
                    Vector3f col(0,0,0);

                    for (int sampleIndex = waveStart; sampleIndex < waveEnd; sampleIndex++) {
                        Float u = Float(i + sampler->sample()) / Float(nx);
                        Float v = Float(j + sampler->sample()) / Float(ny);

                        Ray r = camera->get_ray(u, v);
                        col += color(r, *aggregate, 4);
                    }
                    col /= (waveEnd - waveStart);
                    col = Vector3f(
                            std::sqrt(col.x),
                            std::sqrt(col.y),
                            std::sqrt(col.z));
                    std::cout
                        << int(255.99 * col.x) << " "
                        << int(255.99 * col.y) << " "
                        << int(255.99 * col.z) << " "
                        << "\n";
                }
            }

            // Update start and end wave
            waveStart = waveEnd;
            waveEnd = std::min(samplesPerPixel, waveEnd + nextWaveSize);
            nextWaveSize = std::min(2*nextWaveSize, 64);
        }
    }
}
