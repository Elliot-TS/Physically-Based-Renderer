#include "pbrt/integrator.h"
#include "pbrt/aggregates.h"
#include "pbrt/util/parallel.h"

namespace pbrt {
Vector3f color(
    const Ray &r, const Primitive &shape, int bounces
)
{
  auto si = shape.Intersect(r);
  if (si) {
    Ray scattered;
    Vector3f attenuation;
    if (bounces > 0 &&
        si->interaction.material->scatter(
            r, si->interaction, attenuation, scattered
        ))
    {
      return HorizontalProduct(
          attenuation, color(scattered, shape, bounces + 1)
      );
    }
    else
      return Vector3f(0, 0, 0);
  }
  Vector3f unit_direction = Normalize(r.direction);
  Float t = 0.5 * (unit_direction.y + 1.0);
  return (
      (1.0 - t) * Vector3f(0.88, 0.96, 0.96) +
      t * Vector3f(0.45, 0.73, 0.77)
  );
}

/** Integrator **/
std::optional<ShapeIntersection> Integrator::Intersect(
    const Ray &ray, Float tMax
) const
{
  if (aggregate != nullptr) {
    return aggregate->Intersect(ray, tMax);
  }
  else {
    return {};
  }
}
bool Integrator::IntersectP(const Ray &ray, Float tMax) const
{
  if (aggregate != nullptr) {
    return aggregate->IntersectP(ray, tMax);
  }
  else {
    return false;
  }
}

/** ImageTileIntegrator **/
void ImageTileIntegrator::Render()
{
  bool quit = false;

  Sampler *sampler =
      samplerPrototype;  // We're going to use a ScratchBuffer
                         // later, that's why we're making
                         // this alias

  int spp = 128;  // Samples per pixel should be part of sampler

  int wid = camera->film->display->width;
  int hei = camera->film->display->height;

  // Divide the rendering among a number of threads
  ParallelFor2D(
      wid, hei,
      [&](const unsigned int startIndex,
          const unsigned int endIndex,
          const uint8_t threadID) {
        // Render the image in waves
        int waveSampleStart = 0;
        int waveSampleEnd = 1;
        while (waveSampleStart < spp &&
               camera->film->display->isOpen)
        {
          // Each thread goes through its assigned pixels
          unsigned int x = startIndex % wid;
          unsigned int y = startIndex / wid;
          for (int i = startIndex; i < endIndex; i++) {
            x++;
            if (x >= wid) {
              x = 0;
              y++;
            }

            Vector3f col(0, 0, 0);

            // Compute the samples for each pixel for each
            // wave
            unsigned int sample = waveSampleStart;
            for (; sample < std::min(waveSampleEnd, spp);
                 ++sample) {
              Float u =
                  Float(x + sampler->sample()) / Float(wid);
              Float v =
                  Float(y + sampler->sample()) / Float(hei);

              Ray r = camera->get_ray(u, v);
              col += color(r, *aggregate, 6);
            }
            col /= sample - waveSampleStart;

            col = Vector3f(
                std::sqrt(col.x),
                std::sqrt(col.y),
                std::sqrt(col.z)
            );

            camera->film->AddSample(
                col, x, y, waveSampleStart,
                sample - waveSampleStart
            );

            if (!camera->film->display->isOpen) break;
          }  // End pixel for loop
          if (!camera->film->display->isOpen) break;
          waveSampleStart = waveSampleEnd;
          waveSampleEnd = std::min(waveSampleEnd * 2, spp);
        }  // End waves while loop
      }  // End ParallelFor2D
  );
}
/*void ImageTileIntegrator::Render() {
// TODO: Declare common variables for rendering image in tiles
// TODO: Multi-threading
int nx = camera->film->display->width;
int ny = camera->film->display->height;
Uint32 *imageData = camera->film->display->GetImageData();

Sampler *sampler = samplerPrototype; // We're going to use a
ScratchBuffer later, that's why we're making this alias

// Render image in waves
int samplesPerPixel = 256;
int waveStart = 0, waveEnd = samplesPerPixel, nextWaveSize =
1; while (waveStart < samplesPerPixel)
{
// Render current wave's image tiles in parallel
// TODO: Use mutli-threading via ParallelFor2D() function. For
now, assume there is just one tile.

for (int y = 0; y < ny; ++y) {
for (int x = 0; x < nx; ++x) {
Vector3f col(0,0,0);

for (int sampleIndex = waveStart; sampleIndex < waveEnd;
sampleIndex++) { Float u = Float(x + sampler->sample()) /
Float(nx); Float v = Float(y + sampler->sample()) / Float(ny);

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
}  // namespace pbrt
