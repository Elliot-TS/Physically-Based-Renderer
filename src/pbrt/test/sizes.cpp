#include <iostream>
#include <string>
#include "pbrt/math/vecmath.h"
#include "pbrt/shapes.h"
#include "pbrt/primitive.h"
using namespace pbrt;

#define PrintSize(x,n) std::cout << #x << std::string(1+n, '\t') << sizeof(x) << std::endl

void PrintLine() {
    std::cout << "----------------------------" << std::endl;
}

void MathSizes() {
    PrintLine();
    PrintSize(Float, 2);
    PrintSize(Vector3f, 1);
    PrintSize(Point3f, 2);
}

void PrimitiveSizes() {
    PrintLine();
    PrintSize(Shape, 2);
    PrintSize(Sphere, 2);
    PrintSize(Primitive, 1);
    PrintSize(GeometricPrimitive, 0);
    PrintSize(SimpleAggregate, 1); 
    PrintSize(Material, 1);
}

int main(int argc, char* argv[]) {
    std::cout << "Class name\t\tSize" << std::endl;
    MathSizes();
    PrimitiveSizes();
    return 0;
}
