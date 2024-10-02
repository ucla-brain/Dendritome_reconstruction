//
// Created by muyezhu on 3/1/19.
//
#include <cstring>
#include <string>
#include <regex>
#include <omp.h>
#include "common/mcp3d_common.hpp"
#include "mcp3d_tiff_utils.hpp"

using namespace std;

mcp3d::TiffDirectoryInfo::TiffDirectoryInfo(const string& tiff_path_, int n): tiff_path(tiff_path_)
{
    TIFF* tiff = mcp3d::OpenTiffPath(tiff_path_);
    GetTiffDirectoryInfo(tiff, n);
    TIFFClose(tiff);
}

mcp3d::TiffDirectoryInfo::TiffDirectoryInfo(TIFF* tiff, int n)
{
    MCP3D_ASSERT(tiff)
    tiff_path = string(TIFFFileName(tiff));
    GetTiffDirectoryInfo(tiff, n);
}

void mcp3d::TiffDirectoryInfo::GetTiffDirectoryInfo(TIFF *tiff, int n)
{
    if (!tiff)
        MCP3D_RUNTIME_ERROR("bad pointer to tif image");
    tdir_t current_directory = TIFFCurrentDirectory(tiff);

    // count directories
    n_directories = TIFFNumberOfDirectories(tiff);
    directory_number = n;
    // tiff info of directory n
    MCP3D_ASSERT(n >= 0 && n < n_directories)
    if (!TIFFSetDirectory(tiff, (uint16_t)n))
        MCP3D_RUNTIME_ERROR("error setting tiff directory to" + to_string(n))
    // get tags common to both striped and tiled tiff
    mcp3d::GetCommonTags(tiff, &image_height, &image_width, &bits_per_sample, &samples_per_pixel, &planar_config, description);
    if (TIFFIsTiled(tiff))
    {
        is_tiled = true;
        TIFFGetField(tiff, TIFFTAG_TILEWIDTH, &tile_width);
        TIFFGetField(tiff, TIFFTAG_TILELENGTH, &tile_height);
        rows_in_strip = 0;
    }
    else
    {
        is_tiled = false;
        TIFFGetField(tiff, TIFFTAG_ROWSPERSTRIP, &rows_in_strip);
        tile_height = 0;
        tile_width = 0;
    }
    // guarding against incorrectly written samples_per_pixel fields
    if (samples_per_pixel != 1 && samples_per_pixel != 3 && samples_per_pixel != 4)
    {
        cout << "warning: invalid samples_per_pixel = " << samples_per_pixel << " encountered. calculating from tile / strip size. ";
        if (is_tiled)
        {
            long n_bytes = TIFFTileRowSize(tiff);
            samples_per_pixel = (short)(n_bytes / tile_width / (bits_per_sample / 8));
        }
        else
        {
            long n_bytes = TIFFStripSize(tiff);
            samples_per_pixel = (short)(n_bytes / rows_in_strip / image_width / (bits_per_sample / 8));
        }
        cout << "samples_per_pixel = " << samples_per_pixel << endl;
    }
    // restore directory to value at function entry
    TIFFSetDirectory(tiff, current_directory);
}

void mcp3d::TiffDirectoryInfo::ShowInfo()
{
    cout << "current directory number: " + to_string(directory_number) + "\n" +
            "image width: " + to_string(image_width) + "\n" +
            "image length: " + to_string(image_height) + "\n" +
            "tile width: " + to_string(tile_width) + "\n" +
            "tile width: " + to_string(tile_height) + "\n" +
            "rows in strip: " + to_string(rows_in_strip) + "\n" +
            "bits per sample:" + to_string(bits_per_sample) + "\n" +
            "samples per pixel: " + to_string(samples_per_pixel) + "\n" +
            "planar config: " + to_string(planar_config) + "\n" +
            "directory number: " + to_string(n_directories) + "\n" +
            "is tiled: " + to_string(is_tiled) + "\n" +
            "image description: " + description << endl;
}

bool mcp3d::TiffDirectoryInfo::operator==(const mcp3d::TiffDirectoryInfo& other) const
{
    return image_width == other.image_width &&
           image_height == other.image_height &&
           tile_height == other.tile_height &&
           tile_width == other.tile_width &&
           rows_in_strip == other.rows_in_strip &&
           bits_per_sample == other.bits_per_sample &&
           samples_per_pixel == other.samples_per_pixel &&
           planar_config == other.planar_config &&
           n_directories == other.n_directories &&
           directory_number == other.directory_number &&
           planar_config == other.planar_config &&
           description == other.description;
}

void mcp3d::TiffDirectoryInfo::set_directory_number(int dir_number)
{
    MCP3D_ASSERT(dir_number >= 0)
    directory_number = dir_number;
}

::TIFF* mcp3d::OpenTiffPath(const std::string &tiff_path)
{
    if (!mcp3d::IsFile(tiff_path))
        MCP3D_OS_ERROR(tiff_path + " is not a file");
    TIFF* tiff = TIFFOpen(tiff_path.c_str(), "r");
    if (!tiff)
        MCP3D_RUNTIME_ERROR(tiff_path + " can not be opened as tiff image")
    return tiff;
}

void mcp3d::SetCommonTags(::TIFF *tiff)
{
    MCP3D_ASSERT(tiff)
    TIFFSetField(tiff, TIFFTAG_PLANARCONFIG, 1);
    TIFFSetField(tiff, TIFFTAG_ORIENTATION, 1);
    TIFFSetField(tiff, TIFFTAG_PHOTOMETRIC, 1);
    TIFFSetField(tiff, TIFFTAG_COMPRESSION, COMPRESSION_LZW);
}

void mcp3d::SetTiledTiffTags(::TIFF *tiff, uint32_t image_height, uint32_t image_width, uint32_t tile_height, uint32_t tile_width,
                             short bits_per_sample, short samples_per_pixel)
{
    MCP3D_ASSERT(tiff)
    TIFFSetField(tiff, TIFFTAG_IMAGEWIDTH, image_width);
    TIFFSetField(tiff, TIFFTAG_IMAGELENGTH, image_height);
    TIFFSetField(tiff, TIFFTAG_TILEWIDTH, tile_width);
    TIFFSetField(tiff, TIFFTAG_TILELENGTH, tile_height);
    TIFFSetField(tiff, TIFFTAG_BITSPERSAMPLE, bits_per_sample);
    TIFFSetField(tiff, TIFFTAG_SAMPLESPERPIXEL, samples_per_pixel);
    mcp3d::SetCommonTags(tiff);
}

void mcp3d::SetStripTiffTags(::TIFF *tiff, uint32_t image_height, uint32_t image_width, uint32_t rows_in_strip, short bits_per_sample, short samples_per_pixel)
{
    MCP3D_ASSERT(tiff)
    TIFFSetField(tiff, TIFFTAG_IMAGEWIDTH, image_width);
    TIFFSetField(tiff, TIFFTAG_IMAGELENGTH, image_height);
    TIFFSetField(tiff, TIFFTAG_ROWSPERSTRIP, rows_in_strip);
    TIFFSetField(tiff, TIFFTAG_BITSPERSAMPLE, bits_per_sample);
    TIFFSetField(tiff, TIFFTAG_SAMPLESPERPIXEL, samples_per_pixel);
    mcp3d::SetCommonTags(tiff);
}

void mcp3d::SetOmeTiffTags(::TIFF *tiff, uint32_t image_height, uint32_t image_width, uint32_t tile_height, uint32_t tile_width,
                           int volume_zdim, int volume_ydim, int volume_xdim, short bits_per_sample, short samples_per_pixel)
{
    mcp3d::SetTiledTiffTags(tiff, image_height, image_width, tile_height, tile_width, bits_per_sample, samples_per_pixel);
    TIFFSetField(tiff, TIFFTAG_IMAGEDESCRIPTION, mcp3d::OmeTiffDescription(volume_zdim, volume_ydim, volume_xdim));
}

void mcp3d::SetTiffHyperStackTags(::TIFF *tiff, int n_channels, int n_planes, int n_times, short samples_per_pixel,
                                  short bits_per_sample, float resolution, float min_val, float max_val)
{
    MCP3D_ASSERT(tiff);
    if (samples_per_pixel != 1)
        MCP3D_DOMAIN_ERROR("tiff hyper stack expect samples per pixel = 1")
    MCP3D_ASSERT(n_channels > 0 && n_planes > 0 && n_times > 0)
    string metadata = "ImageJ=1.50e\n";
    int n_imgs = n_channels * n_planes * n_times;
    metadata.append("images=").append(to_string(n_imgs)).append("\n");
    metadata.append("channels=").append(to_string(n_channels)).append("\n");
    metadata.append("slices=").append(to_string(n_planes)).append("\n");
    metadata.append("frames=").append(to_string(n_times)).append("\n");
    metadata.append("hyperstack=true\nmode=color\nunit=micron\nloop=false");
    if (resolution > 0)
        metadata.append("spacing=").append(to_string(resolution)).append("\n");
    if (max_val < 0)
    {
        if (bits_per_sample == 8)
            max_val = float(UINT8_MAX);
        else
            max_val = float(UINT16_MAX);
    }
    metadata.append("min=").append(to_string(min_val)).append("\n");
    metadata.append("max=").append(to_string(max_val)).append("\n");
    TIFFSetField(tiff, TIFFTAG_IMAGEDESCRIPTION, metadata);
}

string mcp3d::GetImageDescriptionTag(::TIFF *tiff)
{
    MCP3D_ASSERT(tiff)
    // get image description string
    unique_ptr<char[]> description_chars;
    TIFFGetField(tiff, TIFFTAG_IMAGEDESCRIPTION, description_chars.get());
    return description_chars ? string(description_chars.get()) : string{};
}

void mcp3d::GetCommonTags(::TIFF *tiff, uint32_t *image_height, uint32_t *image_width, short *bits_per_sample,
                          short *samples_per_pixel, short *planar_config, string& description)
{
    MCP3D_ASSERT(tiff)
    TIFFGetField(tiff, TIFFTAG_IMAGEWIDTH, image_width);
    TIFFGetField(tiff, TIFFTAG_IMAGELENGTH, image_height);
    TIFFGetField(tiff, TIFFTAG_BITSPERSAMPLE, bits_per_sample);
    TIFFGetField(tiff, TIFFTAG_SAMPLESPERPIXEL, samples_per_pixel);
    TIFFGetField(tiff, TIFFTAG_PLANARCONFIG, &planar_config);
    description = mcp3d::GetImageDescriptionTag(tiff);
}

bool mcp3d::GetOmeTiffVolumeHeightWidth(const string &tiff_path, int& zdim, int& ydim, int& xdim)
{
    ::TIFF *tiff = mcp3d::OpenTiffPath(tiff_path);
    string image_description = mcp3d::GetImageDescriptionTag(tiff);
    TIFFClose(tiff);
    regex ome_pattern{mcp3d::OmeTiffDescriptionPattern()};
    smatch m;
    bool found = regex_search(image_description, m, ome_pattern);
    zdim = found ? stoi(m.str(1)) : 0;
    ydim = found ? stoi(m.str(2)) : 0;
    xdim = found ? stoi(m.str(3)) : 0;
    return found;
}

mcp3d::TiffDirectoryInfo mcp3d::VerifyTiffSequence(const string &channel_pyr_dir, const unordered_map<string, vector<string>> &slice_image_names)
{
    mcp3d::TiffDirectoryInfo tiff_info_ref;
    if (slice_image_names.empty())
        return tiff_info_ref;
    // obtain tiff info from one image
    vector<string> slice_names;
    for (const auto& item: slice_image_names)
    {
        slice_names.push_back(item.first);
        if (tiff_info_ref.image_height == 0 && !item.second.empty())
            tiff_info_ref = mcp3d::TiffDirectoryInfo(mcp3d::JoinPath(channel_pyr_dir, item.first, item.second[0]));
    }
    if (tiff_info_ref.image_height == 0)
        return tiff_info_ref;

    // check all tiff images under channel_pyr_dir have same height, width, pixel
    // type and sample number
    bool is_mpi = mcp3d::MPIInitialized();
    int n_threads = is_mpi ? 1 : min(mcp3d::DefaultNumThreads(), (int)slice_image_names.size());
    int n_slices_per_thread = (int)slice_image_names.size() / n_threads;
    mcp3d::MultiThreadExceptions me;
    #pragma omp parallel num_threads(n_threads)
    {
        me.RunAndCaptureException([] {CHECK_PARALLEL_MODEL});
        if (me.HasCapturedException())
        {
        #pragma omp cancel parallel
        }
        me.RunAndCaptureException([&]
        {
             int thread_id = omp_get_thread_num();
             int slice_begin = thread_id * n_slices_per_thread,
                 slice_end = min(slice_begin + n_slices_per_thread, (int)slice_names.size());
             for (int slice_id = slice_begin; slice_id < slice_end; ++slice_id)
             {
                 string slice_name = slice_names[slice_id];
                 const vector<string>& tiff_names = slice_image_names.at(slice_name);
                 for (const auto& tiff_name:tiff_names)
                 {
                     string tiff_path(mcp3d::JoinPath(channel_pyr_dir, slice_name, tiff_name));
                     mcp3d::TiffDirectoryInfo tiff_info(tiff_path);
                     if (tiff_info != tiff_info_ref)
                         MCP3D_RUNTIME_ERROR(tiff_path + " and " + tiff_info_ref.tiff_path + " have different dimensions, pixel data type or samples per pixel")
                 }
             }
        });
        if (me.HasCapturedException())
        {
        #pragma omp cancel parallel
        }
    }
    if (me.HasCapturedException())
        MCP3D_RETHROW(me.e_ptr())
    return tiff_info_ref;
}

mcp3d::VoxelType mcp3d::TiffBitsPerSampleToVoxelType(int bits_per_sample)
{
    if (bits_per_sample == 8)
        return mcp3d::VoxelType::M8U;
    else if (bits_per_sample == 16)
        return mcp3d::VoxelType::M16U;
    else if (bits_per_sample == 32)
        return mcp3d::VoxelType::M32U;
    else if (bits_per_sample == 64)
        return mcp3d::VoxelType::M64U;
    else
        return mcp3d::VoxelType::UNKNOWN;
}
