//
// Created by mzhu on 12/4/18.
//
#include "image/mcp3d_image.hpp"
#include "image/mcp3d_tiff_format.hpp"

using namespace std;

int main(int argc, char** argv)
{
    string tif_dir("/ifs/loni/faculty/dong/mcp/muye/projects/mcp3d/benchmark/profile/tiff_format");
    mcp3d::MImage img(true);
    img.ReadImageInfo(tif_dir);
    mcp3d::MTiffFormat tiff_format {};
    tiff_format.WriteImagePyramidOMP(img, 0, false);
}
