#include "pbrt/integrator.h"
#include "pbrt/util/parallel.h"

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
        bool quit = false;

        Sampler *sampler = samplerPrototype; // We're going to use a ScratchBuffer later, that's why we're making this alias

        int wid = camera->film->display->width;
        int hei = camera->film->display->height;

        // TODO: Render in waves
        ParallelFor2D(wid, hei, 
                [&](const unsigned int startIndex,
                    const unsigned int endIndex,
                    const uint8_t threadID)
                {
                    int samples = 1;
                    while (samples <= 64 && camera->film->display->isOpen) {
                        unsigned int x = startIndex % wid;
                        unsigned int y = startIndex / wid;
                        //while (camera->film->display->isOpen) {
                        for (int i = startIndex; i < endIndex; i++) {
                            x++;
                            if (x >= wid) { x = 0; y++; }

                            Vector3f col(0,0,0);

                            for (unsigned int sample = 0; sample < samples; sample++) {
                                Float u = Float(x + sampler->sample()) / Float(wid);
                                Float v = Float(y + sampler->sample()) / Float(hei);

                                Ray r = camera->get_ray(u, v);
                                col += color(r, *aggregate, 4);
                            }
                            col /= samples;

                            col = Vector3f(
                                    std::sqrt(col.x),
                                    std::sqrt(col.y),
                                    std::sqrt(col.z));
                            camera->film->AddSample(col, x,y, samples-1, samples);

                            if (!camera->film->display->isOpen) break;
                        }
                        if (!camera->film->display->isOpen) break;
                        samples *= 2;
                    }
                });
                }
        /*void ImageTileIntegrator::Render() {
        // TODO: Declare common variables for rendering image in tiles
        // TODO: Multi-threading
        int nx = camera->film->display->width;
        int ny = camera->film->display->height;
        Uint32 *imageData = camera->film->display->GetImageData();

        Sampler *sampler = samplerPrototype; // We're going to use a ScratchBuffer later, that's why we're making this alias

        // Render image in waves
        int samplesPerPixel = 256;
        int waveStart = 0, waveEnd = samplesPerPixel, nextWaveSize = 1;
        while (waveStart < samplesPerPixel)
        {
        // Render current wave's image tiles in parallel
        // TODO: Use mutli-threading via ParallelFor2D() function.  For now, assume there is just one tile.

        for (int y = 0; y < ny; ++y) {
        for (int x = 0; x < nx; ++x) {
        Vector3f col(0,0,0);

        for (int sampleIndex = waveStart; sampleIndex < waveEnd; sampleIndex++) {
        Float u = Float(x + sampler->sample()) / Float(nx);
        Float v = Float(y + sampler->sample()) / Float(ny);

        Ray r = camera->get_ray(u, v);
        col += color(r, *aggregate, 4);
        }
        col /= (waveEnd - waveStart);
        col = Vector3f(
        std::sqrt(col.x),
        std::sqrt(col.y),
        std::sqrt(col.z));
        imageData[y*nx + x] = 
        (Uint8(255.99*col.x) << 24) | 
        (Uint8(255.99*col.y) << 16) | 
        (Uint8(255.99*col.z) << 8) | 
        255;
        }
        }

        // Update start and end wave
        waveStart = waveEnd;
        waveEnd = std::min(samplesPerPixel, waveEnd + nextWaveSize);
        nextWaveSize = std::min(2*nextWaveSize, 64);
        }
        }*/
    }
