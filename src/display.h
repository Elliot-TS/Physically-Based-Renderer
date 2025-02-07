#pragma once
#include "pbrt/math/vecmath.h"
#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
#include <string>

class Image {
private:
   const int xRes;
   const int yRes;
   cv::Mat imageData;

   static std::string windowName;

public:
    struct Color {
        uchar red;
        uchar green;
        uchar blue;
    };

    Image(const int xRes, const int yRes): xRes(xRes), yRes(yRes)
    { imageData = cv::Mat(yRes, xRes, CV_8UC3, cv::Scalar(0,0,0)); }

    // ux varies from 0 (left) to 1 (right)
    // uy varies from 0 (bottom) to 1 (top)
    // Ranges are inclusive: [0,1]
    void createImage ( std::function<Image::Color(const float ux, const float uy)>);

    void Display() {
        cv::namedWindow(windowName, cv::WINDOW_GUI_NORMAL);
        cv::imshow(windowName, imageData);
    }

};
