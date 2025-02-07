#include <iostream>
#include "display.h"
#include <math.h>
//#include <chrono>

std::string Image::windowName = "Image";

int main(int argc, char* argv[]) {
    //using std::chrono::high_resolution_clock;
    //using std::chrono::duration_cast;
    //using std::chrono::duration;
    //using std::chrono::milliseconds;

    const int nx = 3840;
    const int ny = 2160;
    Image img(nx, ny);
    std::cout << "P3\n" << nx << " " << ny << "\n255\n";

    long loopCount = 0;
    long maxLoops = 1000;
    //auto t1 = high_resolution_clock::now();
    while (loopCount == 0) {
        loopCount++;
        img.createImage([loopCount](const float ux, const float uy) -> Image::Color
        {
            return {
                .red=uchar(ux*255),
                .green=uchar(uy*255),
                .blue=uchar(255*std::sin(loopCount/100.0))
            };
        });

        img.Display();
    }
    //auto t2 = high_resolution_clock::now();
    //auto ms_int = duration_cast<milliseconds>(t2-t1);
    //duration<double, std::milli> ms_double = t2 - t1;
    //std::cout << ms_int.count() << "ms\n";
    //std::cout << ms_double.count() << "ms\n";
    
    return 0;
}
