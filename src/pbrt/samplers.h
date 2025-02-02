#pragma once
#include <random>
#include "util/float.h"
#include "util/vecmath.h"

namespace pbrt {
    class Sampler{
        private:
            std::mt19937 generator; 

        public:
            virtual inline Float sample() = 0;
            virtual inline Vector3f s_Vector3f() = 0;
            virtual inline Point3f s_uSphere() = 0;
    };

    class UniformSampler : public Sampler{
        private:
            std::uniform_real_distribution<Float> distribution;
            std::mt19937 generator; 

        public:
            UniformSampler() {
                // obtain a random number from hardware
                std::random_device rd; 
                // seed the generator
                std::mt19937 generator(rd()); 
                // define the range
                std::uniform_real_distribution<Float> distribution(0.0, 1.0); 
            }

            inline Float sample() {
                return distribution(generator); 
            }
            inline Vector3f s_Vector3f() {
                return Vector3f(sample(), sample(), sample());
            }
            inline Point3f s_uSphere()
            {
                Vector3f p;
                do {
                    p = 2.0*s_Vector3f() - Vector3f(1,1,1);
                } while (LengthSquared(p) >= 1.0);

                return Point3f(p);
            }
    };
}
