//
// Created by muyezhu on 3/1/19.
//

#ifndef MCP3D_MCP3D_TIFF_UTILS_HPP
#define MCP3D_MCP3D_TIFF_UTILS_HPP

#include <unordered_map>
#include <tiff.h>
#include <tiffio.h>
#include "mcp3d_voxel_types.hpp"

namespace mcp3d
{

struct TiffDirectoryInfo
{
    TiffDirectoryInfo() : image_width(0), image_height(0), tile_width(0), tile_height(0), rows_in_strip(0),
                          bits_per_sample(0), samples_per_pixel(0), planar_config(-1), n_directories(0),
                          directory_number(-1), is_tiled(true), tiff_path(std::string{}) {};
    /// obtain MCPTIFFInfo fields.
    /// if multiple directory exists for the tif image,
    /// n_directory will be updated accordingly,
    /// but other fields will use values of first directory
    /// tiff file handle created in the constructor will be closed
    explicit TiffDirectoryInfo(const std::string& tiff_path_, int n = 0);

    /// does not close the tiff handle
    explicit TiffDirectoryInfo(::TIFF* tiff, int n = 0);

    void ShowInfo();

    /// records the current directory of tiff. obtain requested information, and
    /// sets current directory back to original
    void GetTiffDirectoryInfo(TIFF *tiff, int n = 0);

    bool empty() const  { return image_width == 0; }

    // return true if image_height, image_width, bits_per_sample,
    // samples_per_pixel, n_directories and planar_config are identical
    // between two TiffDirectoryInfo structs
    bool operator==(const TiffDirectoryInfo &other) const;

    bool operator!=(const TiffDirectoryInfo &rhs) const
    { return !((*this) == rhs); }

    /// note: only sets directory_number, other fields under new directory number
    /// will not be examined. will do so for empty instances as well. intended
    /// use case is with
    /// ReadTiffDirectoryData(TIFF* tiff, const TiffDirectoryInfo& tiff_dir_info,
    /// int n, std::unique_ptr<uint8_t[]>& data), when caller is certain all
    /// directory has same fields other than directory number
    void set_directory_number(int dir_number);

    uint32_t image_width, image_height;
    uint32_t tile_width, tile_height;
    uint32_t rows_in_strip;
    short bits_per_sample, samples_per_pixel;
    short planar_config;
    short n_directories;
    int directory_number;
    bool is_tiled;
    std::string tiff_path, description;
};

::TIFF* OpenTiffPath(const std::string& tiff_path);

/// set TIFFTAG_PLANARCONFIG = 1, TIFFTAG_ORIENTATION = 1,
/// TIFFTAG_PHOTOMETRIC = 1. these tags must be set for opencv to correctly
/// read an image
void SetCommonTags(::TIFF *tiff);

void SetTiledTiffTags(::TIFF *tiff, uint32_t image_height, uint32_t image_width, uint32_t tile_height, uint32_t tile_width,
                      short bits_per_sample, short samples_per_pixel);

void SetStripTiffTags(::TIFF *tiff, uint32_t image_height, uint32_t image_width, uint32_t rows_in_strip, short bits_per_sample, short samples_per_pixel);

/// mcp3d ome tiff image description
/// volume height: xxx
// volume width: xxx
//// volume depth: xxx
inline std::string OmeTiffDescription(int zdim, int ydim, int xdim)
{ return "mcp3d ome tiff image description\n"
         "volume depth: " + std::to_string(zdim) + "\nvolume height: " + std::to_string(ydim) + "\nvolume width: " + std::to_string(xdim); }

inline std::string OmeTiffDescriptionPattern()
{ return "mcp3d ome tiff image description\nvolume depth: ([0-9]+)\nvolume height: ([0-9]+)\nvolume width: ([0-9]+)\n"; }

/// volume_dim is the dimension of the volume composed of small ometiff images.
/// the ometiff images will be written as tiled. tile_height is the height of the tiles
void SetOmeTiffTags(::TIFF *tiff, uint32_t image_height, uint32_t image_width, uint32_t tile_height, uint32_t tile_width,
                    int volume_zdim, int volume_ydim, int volume_xdim, short bits_per_sample, short samples_per_pixel);

/*
 * see ImageJ hyperstack plugin for some details:
 * https://imagej.nih.gov/ij/developer/api/ij/plugin/HyperStackConverter.html
 * most relevant: The default "xyczt" order is used if 'order' is null
 * DimensionOrder attribute that specifies the rasterization order of the
 * image planes. For example, XYZTC means that there is a series of image
 * planes with the Z axis varying fastest, followed by T, followed by C
 * (e.g. if a XYZTC dataset contains two focal planes, three time points and
 * two channels, the order would be: Z0-T0-C0, Z1-T0-C0, Z0-T1-C0, Z1-T1-C0,
 * Z0-T2-C0, Z1-T2-C0, Z0-T0-C1, Z1-T0-C1, Z0-T1-C1, Z1-T1-C1, Z0-T2-C1,
 * Z1-T2-C1).
 * example hyperstack tiff_info output:
 * TIFF Directory at offset 0xe233bb25 (3795041061)
 * Subfile Type: (0 = 0x0)
 * Image Width: 1400 Image Length: 1400
 * Resolution: 3.55544, 3.55544 (unitless)
 * Bits/Sample: 16
 * Compression Scheme: None
 * Photometric Interpretation: min-is-black
 * Samples/Pixel: 1
 * Rows/Strip: 1400
 * Planar Configuration: single image plane
 * ImageDescription: ImageJ=1.50e
 * images=968  // this is channels * slices * time points
 * channels=4
 * slices=242  // this is z level
 * hyperstack=true
 * mode=color
 * unit=micron
 * spacing=3.0
 * loop=false
 * min=0.0
 * max=65535.0
 * The whole ImageDescription string is the OME-XML portion of OME-TIF
*/
void SetTiffHyperStackTags(::TIFF *tiff, int n_channels, int n_planes, int n_times = 1, short samples_per_pixel = 1,
                           short bits_per_sample = 8, float resolution = -1.0f, float min_val = 0.0f, float max_val = -1.0f);

std::string GetImageDescriptionTag(::TIFF *tiff);

void GetCommonTags(::TIFF *tiff, uint32_t* image_height, uint32_t* image_width, short* bits_per_sample,
                   short* samples_per_pixel, short* planar_config, std::string& description);

/// return true if OmeTiffDescriptionPattern() found in field TIFFTAG_IMAGEDESCRIPTION, parse volume dims
/// otherwise set volume dims to 0 and return false
bool GetOmeTiffVolumeHeightWidth(const std::string& tiff_path, int& zdim, int& ydim, int& xdim);

/// slice_image_names is in the same organization as MChannelPyrInfo
/// slice_image_names_ member. assert that all images have identical TiffDirectoryInfo
/// and returns the TiffDirectoryInfo
TiffDirectoryInfo VerifyTiffSequence(const std::string& channel_pyr_dir, const std::unordered_map<std::string, std::vector<std::string>>& slice_image_names);

mcp3d::VoxelType TiffBitsPerSampleToVoxelType(int bits_per_sample);

}


#endif //MCP3D_MCP3D_TIFF_UTILS_HPP
