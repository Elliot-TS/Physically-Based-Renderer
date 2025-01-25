#include <iostream>
#include <string>
#include "pbrt/util/vecmath.h"

using namespace pbrt;

int main(int arcg, char *argv[]) {
    Vector2<bool> t1;
    Vector2<int> t2(5,3);
    Vector2<float> t3(5.32f, -43.f);
    Vector2<std::string> t4("Hello", "World");

    std::cout << t2 << std::endl;
    t2 += t3;
    std::cout << t2 << std::endl;
    std::cout << Abs(t2) << std::endl;


    // Ray tracing in a weekend
    //int nx = 200;
    //int ny = 100;
    //std::cout << "P3\n" << nx << " " << ny << "\n255\n";

    //for (int j = ny-1; j >= 0; j--) {
        //for (int i = 0; i < nx; i++) {
            //Vector3<float> col(
                //float(i) / float(nx),
                //float(j) / float(ny),
                //0.2
            //);
            //int ir = int(255.99*col.x);
            //int ig = int(255.99*col.y);
            //int ib = int(255.99*col.z);
            //std::cout << ir << " " << ig << " " << ib << std::endl;
        //}
    //}
        


    // Convert command-line arguments to vector of strings TODO
    // Declare variables for parsed command line TODO
    // Process command-line arguments TODO
    // Initialize pbrt TODO
    // Parse provided scene description files TODO
    // Render the scene TODO
    // Clean up after rendering the scene TODO
}
