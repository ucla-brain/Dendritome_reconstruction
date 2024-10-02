//
// Created by muyezhu on 3/2/18.
//

#ifndef MCP3D_MCP3D_IMAGE_CONSTANTS_HPP
#define MCP3D_MCP3D_IMAGE_CONSTANTS_HPP

namespace mcp3d
{
/// tile width and height in tiffs written as tiled
static const int TIFFTILE_XDIM = 128;
static const int TIFFTILE_YDIM = 128;
/// width and height of ome.tif
static const int TIFFCHUNK_XDIM = 2048;
static const int TIFFCHUNK_YDIM = 2048;
static const int TIFFCHUNK_ZDIM = 100;
static const int ZSCALE_NONE = -1;
/// number of digits in image volume dimensions
static const int VOLUME_DIM_WIDTH = 6;
/// number of digits in ome.tif files chunk dimensions
static const int OMETIFF_DIM_WIDTH = 6;

}

#endif //MCP3D_MCP3D_IMAGE_CONSTANTS_HPP
