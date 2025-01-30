#pragma once
#include <vector>
class Integrator {
    public:
        // Methods
        Integrator(Primitive aggregate, std::vector<Light> lights)
            : aggregate(aggregate), lights(lights)
        {
            // TODO: Integrator Constructor Implementation
        }
        // Members TODO
        Primitive aggregate;
        std::vector<Light> lights;
    protected:
        // Methods TODO
};
