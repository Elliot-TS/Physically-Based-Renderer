#include "display.h"

// See https://docs.opencv.org/2.4/doc/tutorials/core/how_to_scan_images/how_to_scan_images.html#performance-difference
void Image::createImage ( std::function<Image::Color(const float ux, const float uy)> setColor) {
    for (int j = yRes-1; j >= 0; j--) {
        for (int i = 0; i < xRes; i++){
            Image::Color c = setColor(i/float(xRes),j/float(yRes));
            std::cout << int(c.red) << " " << int(c.green) << " " << int(c.blue) << "\n";
        }
    }
    //CV_Assert(imageData.depth() == CV_8U);
    
    //int channels = imageData.channels();
    //int nRows = imageData.rows;
    //int nCols = int(imageData.cols) * channels;

    //if (imageData.isContinuous()) {
        //nCols *= nRows;
        //nRows = 1;
    //}

    //int i,j;
    //uchar* p;
    //for (int i = 0; i < nRows; i++) {
        //p = imageData.ptr<uchar>(i);
        //for (int j = 0; j < nCols; j += channels) {
            //float ux = j % (xRes*channels) 
                /// float((xRes-1)*channels);
            //float uy = (i + int(j/(channels*xRes))) 
                /// float(yRes);
            //Image::Color col = setColor(ux, 1-uy);
            //p[j+0] = col.blue;
            //p[j+1] = col.green;
            //p[j+2] = col.red;
        //}
    //}
}
