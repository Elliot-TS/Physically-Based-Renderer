#include <iostream>
#include <math.h>
#include <chrono>
#include "pbrt/util/display.h"


int main(int argc, char* argv[]) {
    // For Benchmarking
    using std::chrono::high_resolution_clock; using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;

    // Create a window
    pbrt::DisplayWindow display(3840, 2160);

    // Create image buffer
    Uint32 * imageData = display.GetImageData();
    for (int y = 0; y < display.height; y++) {
        for (int x = 0; x < display.width; x++) {
            Uint8 r = x * 255 / display.width;
            Uint8 g = y * 255 / display.height;
            Uint8 b = 0;
            Uint8 a = 255;
            imageData[y*display.width + x] = (r << 24) | (g << 16) | (b << 8) | a;
        }
    }

    display.Open();

    Uint8 blue = 0;
    while (display.isOpen) {
        for (int y = 0; y < display.height; y++) {
            for (int x = 0; x < display.width; x++) {
                Uint8 r = x * 255 / display.width;
                Uint8 g = y * 255 / display.height;
                Uint8 b = blue;
                Uint8 a = 255;
                int count = 0;
                for (int i = 0; i < 10000000; i++) {
                    count += i;
                }
                imageData[y*display.width + x] = (r << 24) | (g << 16) | (b << 8) | a;
            }
        }
        blue++;

    }

    return 0;
}
