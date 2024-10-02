//
// Created by mzhu on 4/4/17.
//

#ifndef MCP3D_LARGE_IMAGE_IO_HPP
#define MCP3D_LARGE_IMAGE_IO_HPP

#include <cstdlib>
#include <iostream>
#include <memory>
#include <tiff.h>
#include <tiffio.h>

/// tdata_t is a typedef of void* in tiffio.h

namespace mcp3d
{

void ResizeLargeTiff(const std::string &input_tiff_path, const std::string& resize_tiff_path);

/// will close tiff handle
void ReadTiffDirectoryData(const std::string& tiff_path, int n, std::unique_ptr<uint8_t[]>& data);

/// will not close tiff handle
/// convert data to gray. if data manages a null pointer and stored iamge is grey,
/// allocate memory for data's pointer
void ReadTiffDirectoryData(::TIFF* tiff, int n, std::unique_ptr<uint8_t[]>& data);

/// will not close tiff handle. the TIFF handle is expected to be set at desired directory
/// if data storage is RGB(A), convert to grey
/// note that for performance reason, if pointer managed by data is not null, and
/// stored image is grey, data will be used directly with no size check
void ReadStripedTiff(::TIFF* tiff, std::unique_ptr<uint8_t[]>& data);

/// will not close tiff handle. the TIFF handle is expected to be set at desired directory
/// if data storage is RGB(A), convert to grey
/// note that for performance reason, if pointer managed by data is not null, and
/// stored image is grey, data will be used directly with no size check
void ReadTiledTiff(::TIFF* tiff, std::unique_ptr<uint8_t[]>& data);

/// note: input can be in unspecified state after function call. output and input
/// can be the same instance
void ColorToGray(std::unique_ptr<uint8_t[]> &input, int zdim, int ydim, int xdim,
                 int bytes_per_sample, int sample_per_pixel, std::unique_ptr<uint8_t[]> &output);

/// copy sub image data from tiled tiff image into subimg_buf
/// the subimg's top left pixel has coordinate (subimg_origin_x, subimg_origin_y)
/// the subimg has height subimg_height and witdh subimg_width
/// if portion of sub image is outside of tiff image, black pixels are padded
/// called by large_image_io::large_tiff_to_subimgs internally
/// if not given or 0, subimg_width is same as subimg_height
/// subimg_src_origin_x: x coordinate of top left pixel of subimage
///                      of original tiff image
/// background: color of background. default 'b': background is black
void TransferTileToSubimg(::TIFF *tif, tdata_t subimg_buf, int64_t subimg_height, int64_t subimg_width,
                          int64_t subimg_src_origin_x, int64_t subimg_src_origin_y);

void TransferTileToSubimg(const std::string& tiff_path, tdata_t subimg_buf, int64_t subimg_height, int64_t subimg_width,
                          int64_t subimg_src_origin_x, int64_t subimg_src_origin_y);

void TransferScanlineToSubimg(const std::string& tiff_path, tdata_t subimg_buf, int64_t subimg_height, int64_t subimg_width,
                              int64_t subimg_src_origin_x, int64_t subimg_src_origin_yund);

void TransferScanlineToSubimg(::TIFF *tif, tdata_t subimg_buf, int64_t subimg_height, int64_t subimg_width,
                              int64_t subimg_src_origin_x, int64_t subimg_src_origin_y);

void TransferToSubimg(::TIFF *tif, tdata_t subimg_buf, int64_t subimg_height, int64_t subimg_width,
                      int64_t subimg_src_origin_x, int64_t subimg_src_origin_y);

void TransferToSubimg(const std::string& tif_path, tdata_t subimg_buf, int64_t subimg_height, int64_t subimg_width,
                      int64_t subimg_src_origin_x, int64_t subimg_src_origin_y);

}

#endif //MCP3D_LARGE_IMAGE_IO_HPP
