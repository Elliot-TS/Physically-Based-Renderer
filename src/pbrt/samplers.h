#pragma once
#include <random>
#include "pbrt/math/vecmath.h"

namespace pbrt {
    class Sampler{
        private:
            std::mt19937 generator; 

        public:
            virtual Float sample() = 0;
            virtual Vector3f s_Vector3f() = 0;
            virtual Point3f s_uSphere() = 0;
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

            Float sample() {
                return distribution(generator); 
            }
            Vector3f s_Vector3f() {
                return Vector3f(sample(), sample(), sample());
            }
            Point3f s_uSphere()
            {
                Vector3f p;
                do {
                    p = 2.0*s_Vector3f() - Vector3f(1,1,1);
                } while (LengthSquared(p) >= 1.0);

                return Point3f(p);
            }
    };
}
