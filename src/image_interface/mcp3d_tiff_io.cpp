//
// Created by mzhu on 4/4/17.
//
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include "mcp3d_tiff_utils.hpp"
#include "mcp3d_tiff_io.hpp"

using namespace std;

void mcp3d::ResizeLargeTiff(const string &input_tiff_path, const string& resize_tiff_path)
{
    ::TIFF* tif = mcp3d::OpenTiffPath(input_tiff_path);
    ::TIFF* tiff_resize = TIFFOpen(resize_tiff_path.c_str(), "w");

    if (tif)
    {
        uint32_t image_width, image_length;
        uint32_t tile_width, tile_length;
        short bits_per_pixel;
        uint32_t x, y;
        tdata_t buf, buf_resize;

        TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &image_width);
        TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &image_length);
        TIFFGetField(tif, TIFFTAG_TILEWIDTH, &tile_width);
        TIFFGetField(tif, TIFFTAG_TILELENGTH, &tile_length);
        TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bits_per_pixel);

        uint32_t img_width_resize = (image_width - image_width % tile_width) / 2;
        uint32_t img_length_resize = (image_length - image_length % tile_length) / 2;
        uint32_t tile_width_resize = tile_width / 2;
        uint32_t tile_length_resize = tile_length / 2;
        mcp3d::SetTiledTiffTags(tiff_resize, img_length_resize, img_width_resize, tile_length_resize, tile_width_resize, bits_per_pixel, 1);

        buf = _TIFFmalloc(TIFFTileSize(tif));
        buf_resize = _TIFFmalloc(TIFFTileSize(tif) / 4);
        ttile_t tile = 0;
        tsize_t size = TIFFTileSize(tif) / 4;

        for (y = 0; y < image_length - tile_length; y += tile_length)
            for (x = 0; x < image_width - tile_width; x += tile_width)
            {
                TIFFReadTile(tif, buf, x, y, 0, 1);
                for (uint32_t i = 0; i < tile_length; ++i)
                    for (uint32_t j = 0; j < tile_width; ++j)
                        if (i % 2 == 0 && j % 2 == 0)
                        {
                            ((uchar*)buf_resize)[i / 2 * tile_width_resize + j / 2] = ((uchar*)buf)[i * tile_width + j];
                        }


                if (TIFFWriteRawTile(tiff_resize, tile, buf_resize, size) > 0)
                    ++tile;
            }

        _TIFFfree(buf);
        TIFFClose(tif);
        TIFFClose(tiff_resize);
        _TIFFfree(buf_resize);
    }
}

void mcp3d::ReadTiffDirectoryData(const string& tiff_path, int n, unique_ptr<uint8_t[]>& data)
{
    ::TIFF* tiff = TIFFOpen(tiff_path.c_str(), "r");
    mcp3d::ReadTiffDirectoryData(tiff, n, data);
    TIFFClose(tiff);
}

void mcp3d::ReadTiffDirectoryData(::TIFF *tiff, int n, std::unique_ptr<uint8_t[]> &data)
{
    MCP3D_ASSERT(tiff)
    MCP3D_ASSERT(TIFFSetDirectory(tiff, (uint16_t)n))
    if (TIFFIsTiled(tiff))
        mcp3d::ReadTiledTiff(tiff, data);
    else
        mcp3d::ReadStripedTiff(tiff, data);
}

void mcp3d::ColorToGray(unique_ptr<uint8_t[]> &input, int zdim, int ydim, int xdim,
                        int bytes_per_sample, int sample_per_pixel, unique_ptr<uint8_t[]> &output)
{
    MCP3D_ASSERT(input)
    MCP3D_ASSERT(zdim > 0 && ydim > 0 && xdim > 0)
    MCP3D_ASSERT(bytes_per_sample == 1 || bytes_per_sample == 2 || bytes_per_sample == 4 || bytes_per_sample == 8)
    MCP3D_ASSERT(sample_per_pixel == 1 || sample_per_pixel == 3 || sample_per_pixel == 4)
    if (sample_per_pixel == 1)
    {
        if (input != output)
            output = move(input);
        return;
    }
    // variable temp is used in case output and input are the same
    unique_ptr<uint8_t[]> temp = move(make_unique<uint8 []>((size_t)ydim * (size_t)xdim * (size_t)bytes_per_sample));
    uint8_t* input_ptr = input.get();
    uint8_t* temp_ptr = temp.get();
    int cv_input_type = mcp3d::MakeCVType(bytes_per_sample, sample_per_pixel),
            cv_output_type = mcp3d::MakeCVType(bytes_per_sample, 1);
    long n_input_plane_bytes = (long)ydim * (long)xdim * (long)sample_per_pixel * (long)bytes_per_sample,
            n_output_plane_bytes = (long)ydim * (long)xdim * (long)bytes_per_sample;
    for (int z = 0; z < zdim; ++z)
    {
        cv::Mat m_in(ydim, xdim, cv_input_type, input_ptr + n_input_plane_bytes * z);
        cv::Mat m_out(ydim, xdim, cv_output_type, temp_ptr + n_output_plane_bytes * z);
        int cvt_mode = sample_per_pixel == 3 ? cv::COLOR_RGB2GRAY : cv::COLOR_RGBA2GRAY;
        SERIALIZE_OPENCV_MPI
        #if CV_MAJOR_VERSION < 4
            cv::cvtColor(m_in, m_out, cvt_mode);
        #else
            cv::cvtColor(m_in, m_out, cvt_mode);
        #endif
    }
    output = move(temp);
}

void mcp3d::ReadStripedTiff(::TIFF *tiff, unique_ptr<uint8_t[]> &data)
{
    MCP3D_ASSERT(tiff)
    short bits_per_sample, samples_per_pixel, planar_config;
    TIFFGetField(tiff, TIFFTAG_BITSPERSAMPLE, &bits_per_sample);
    TIFFGetField(tiff, TIFFTAG_SAMPLESPERPIXEL, &samples_per_pixel);
    TIFFGetField(tiff, TIFFTAG_PLANARCONFIG, &planar_config);
    uint32_t image_width, image_height;
    TIFFGetField(tiff, TIFFTAG_IMAGEWIDTH, &image_width);
    TIFFGetField(tiff, TIFFTAG_IMAGELENGTH, &image_height);
    size_t bytes_per_pixel = (size_t)bits_per_sample / 8 * samples_per_pixel;
    size_t n_image_bytes = (size_t)bytes_per_pixel * image_height * image_width;
    if (!data || samples_per_pixel > 1)
        data = make_unique<uint8_t []>(n_image_bytes);
    uint8_t* data_ptr = data.get();
    tstrip_t n_strips = TIFFNumberOfStrips(tiff);
    tsize_t n_strip_bytes = TIFFStripSize(tiff), data_ptr_offset = 0;
    for (tstrip_t i = 0; (tstrip_t)i < n_strips; ++i)
    {
        TIFFReadEncodedStrip(tiff, i, data_ptr + data_ptr_offset, n_strip_bytes);
        data_ptr_offset += n_strip_bytes;
    }
    if (samples_per_pixel > 1)
        mcp3d::ColorToGray(data, 1, image_height, image_width, bits_per_sample / 8, samples_per_pixel, data);
}

void mcp3d::ReadTiledTiff(::TIFF* tiff, unique_ptr<uint8_t[]> &data)
{
    MCP3D_ASSERT(tiff)
    short bits_per_sample, samples_per_pixel, planar_config;
    TIFFGetField(tiff, TIFFTAG_BITSPERSAMPLE, &bits_per_sample);
    TIFFGetField(tiff, TIFFTAG_SAMPLESPERPIXEL, &samples_per_pixel);
    TIFFGetField(tiff, TIFFTAG_PLANARCONFIG, &planar_config);
    uint32_t image_width, image_height;
    TIFFGetField(tiff, TIFFTAG_IMAGEWIDTH, &image_width);
    TIFFGetField(tiff, TIFFTAG_IMAGELENGTH, &image_height);
    uint32_t tile_width, tile_height;
    TIFFGetField(tiff, TIFFTAG_TILEWIDTH, &tile_width);
    TIFFGetField(tiff, TIFFTAG_TILELENGTH, &tile_height);
    // number of bytes required by image
    size_t bytes_per_pixel = (size_t)bits_per_sample / 8 * samples_per_pixel;
    size_t n_image_bytes = (size_t)bytes_per_pixel * image_height * image_width;
    if (!data || samples_per_pixel > 1)
        data = move(make_unique<uint8_t []>(n_image_bytes));
    uint8_t * data_ptr = data.get();
    tdata_t tile_buffer;
    tile_buffer = _TIFFmalloc(TIFFTileSize(tiff));
    for (uint32_t y = 0; y < image_height; y += tile_height)
        for (uint32_t x = 0; x < image_width; x += tile_width)
        {
            TIFFReadTile(tiff, tile_buffer, x, y, 0, (tsample_t)planar_config);
            // copy data one row at a time
            size_t n_tile_row_bytes = bytes_per_pixel * tile_width;
            size_t data_ptr_offset = bytes_per_pixel * (y * image_width + x);
            size_t tile_buffer_offset = 0;
            for (uint32_t tile_yoffset = 0; tile_yoffset < tile_height; ++ tile_yoffset)
            {
                memcpy(data_ptr + data_ptr_offset, ((uint8_t*)tile_buffer) + tile_buffer_offset, n_tile_row_bytes);
                data_ptr_offset += bytes_per_pixel * image_width;
                tile_buffer_offset += n_tile_row_bytes;
            }
        }
    _TIFFfree(tile_buffer);
    if (samples_per_pixel > 1)
        mcp3d::ColorToGray(data, 1, image_height, image_width, bits_per_sample / 8, samples_per_pixel, data);
}

void ValidateTiffSubImageRequest(const mcp3d::TiffDirectoryInfo& tiff_info, const int64_t subimg_height, const int64_t subimg_width,
                                 const int64_t subimg_origin_x = 0, const int64_t subimg_origin_y = 0)
{
    if (tiff_info.n_directories > 1)
        MCP3D_DOMAIN_ERROR("does not support multi directory tiff image")
    if (tiff_info.planar_config != 1)
        // https://www.awaresystems.be/imaging/tiff/tifftags/planarconfiguration.html
        MCP3D_DOMAIN_ERROR("only supporting chunky planar config")
    if (tiff_info.samples_per_pixel > 3)
        MCP3D_DOMAIN_ERROR("does not support alpha channel")
    if (tiff_info.bits_per_sample != 8 && tiff_info.bits_per_sample != 16 && tiff_info.bits_per_sample != 32)
        MCP3D_DOMAIN_ERROR("only supporting 8 bits, 16 bits or 32 bits per voxel")
    if (subimg_origin_x < 0 || subimg_origin_x >= tiff_info.image_width ||
        subimg_origin_y < 0 || subimg_origin_y >= tiff_info.image_height)
        MCP3D_INVALID_ARGUMENT("invalid subimage starting coordinates")
}

void mcp3d::TransferTileToSubimg(const string& tiff_path, tdata_t subimg_buf, int64_t subimg_height, int64_t subimg_width,
                                 int64_t subimg_src_origin_x, int64_t subimg_src_origin_y)
{
    ::TIFF* tiff = mcp3d::OpenTiffPath(tiff_path);
    mcp3d::TransferTileToSubimg(tiff, subimg_buf, subimg_height, subimg_width, subimg_src_origin_x, subimg_src_origin_y);
    ::TIFFClose(tiff);
}

void mcp3d::TransferTileToSubimg(::TIFF *tif, tdata_t subimg_buf, int64_t subimg_height, int64_t subimg_width,
                                 int64_t subimg_src_origin_x, int64_t subimg_src_origin_y)
{
    MCP3D_ASSERT(subimg_buf)
    mcp3d::TiffDirectoryInfo tiff_info(tif);
    if (! tiff_info.is_tiled)
        MCP3D_DOMAIN_ERROR("wrong function call: tiff image is stripped")
    if (subimg_width == 0)
        subimg_width = subimg_height;
    uint32_t img_height = tiff_info.image_height,
             img_width = tiff_info.image_width,
             tile_height = tiff_info.tile_height,
             tile_width = tiff_info.tile_width;
    short bits_per_sample = tiff_info.bits_per_sample,
          samples_per_pixel = tiff_info.samples_per_pixel,
          planar_config = tiff_info.planar_config;
    ValidateTiffSubImageRequest(tiff_info, subimg_height, subimg_width, subimg_src_origin_x, subimg_src_origin_y);
    tdata_t tile_buf;
    tile_buf = _TIFFmalloc(TIFFTileSize(tif));
    if (!tile_buf)
        MCP3D_BAD_ALLOC("can not allocate memory for image tile")

    // subimg_focus_x, subimg_focus_y:
    // relative to the top left pixel of the original tif image.
    // coordinates of the upper left pixel of rectangular region in the subimg
    // to be copied. does not have to be aligned with
    // top left pixel of the next tile in original tif to be read
    // tile_start_x, tile_start_y:
    // relative to the top left pixel of the next tile to be copied.
    // relative coordinates of subimg_focus point in the tile.
    // needed when subimg dimension do not evenly divide tile dimension
    int64_t tile_start_x, tile_start_y,
            subimg_focus_x, subimg_focus_y;
    subimg_focus_x = subimg_src_origin_x;
    subimg_focus_y = subimg_src_origin_y;
    while (true)
    {
        // fill subimg with current tile.
        // if subimg is partially out of image area, background padding is already
        // done in MImageIO, just make sure to read within bounds
        if (subimg_focus_x < img_width && subimg_focus_y < img_height)
            TIFFReadTile(const_cast<TIFF*>(tif), tile_buf, (uint32_t)subimg_focus_x, (uint32_t)subimg_focus_y, 0, (uint16_t)planar_config);

        // subimg_x_offset, subimg_y_offset:
        // relative to subimg_focus_x and subimg_focus_y,
        // used to correctly index into the rectangular area in the subimg
        // currently being filled
        int64_t subimg_x_offset = subimg_focus_x - subimg_src_origin_x,
                subimg_y_offset = subimg_focus_y - subimg_src_origin_y;
        int64_t subimg_sample_offset =
                samples_per_pixel * (subimg_y_offset * subimg_width + subimg_x_offset);

        tile_start_x = subimg_focus_x % tile_width;
        tile_start_y = subimg_focus_y % tile_height;

        int64_t delta_y = subimg_focus_y % tile_height == 0 ?
                          min((int64)tile_height, (int64)subimg_height - subimg_y_offset) :
                          min((int64)tile_height - subimg_focus_y % tile_height, (int64)subimg_height - subimg_y_offset),
                delta_x = subimg_focus_x % tile_width == 0 ?
                          min((int64)tile_width, (int64)subimg_width - subimg_x_offset) :
                          min((int64)tile_width - subimg_focus_x % tile_width, (int64)subimg_width - subimg_x_offset);
        int64_t tile_end_y = tile_start_y + delta_y, tile_end_x = tile_start_x + delta_x;

        if (subimg_focus_x < img_width && subimg_focus_y < img_height)
        {
            int64_t subimg_sample_index = subimg_sample_offset;
            int64_t tile_sample_index = samples_per_pixel * (tile_width * tile_start_y + tile_start_x);
            size_t n_bytes = (size_t)(tile_end_x - tile_start_x) * (bits_per_sample / 8) * samples_per_pixel;
            for (int64_t i = tile_start_y; i < tile_end_y; ++i)
            {
                if (bits_per_sample == 8)
                    memcpy(((uint8_t*)subimg_buf) + subimg_sample_index, ((uint8_t*)tile_buf) + tile_sample_index, n_bytes);
                else if (bits_per_sample == 16)
                    memcpy(((uint16_t*)subimg_buf) + subimg_sample_index, ((uint16_t*)tile_buf) + tile_sample_index, n_bytes);
                else
                    memcpy(((int32_t*)subimg_buf) + subimg_sample_index, ((int32_t*)tile_buf) + tile_sample_index, n_bytes);
                subimg_sample_index += samples_per_pixel * subimg_width;
                tile_sample_index += samples_per_pixel * tile_width;
            }
        }
        // update subimg_focus_x and subimg_focus_y
        subimg_focus_x += delta_x;
        if (subimg_focus_x - subimg_src_origin_x >= subimg_width)
        {
            subimg_focus_x = subimg_src_origin_x;
            subimg_focus_y += delta_y;
            if (subimg_focus_y - subimg_src_origin_y >= subimg_height)
                break;
        }
    }
    _TIFFfree(tile_buf);
}

void mcp3d::TransferScanlineToSubimg(const string& tiff_path, tdata_t subimg_buf, int64_t subimg_height, int64_t subimg_width,
                                     int64_t subimg_src_origin_x, int64_t subimg_src_origin_y)
{
    ::TIFF* tiff = mcp3d::OpenTiffPath(tiff_path);
    mcp3d::TransferScanlineToSubimg(tiff, subimg_buf, subimg_height, subimg_width, subimg_src_origin_x, subimg_src_origin_y);
    ::TIFFClose(tiff);
}

void mcp3d::TransferScanlineToSubimg(::TIFF *tif, tdata_t subimg_buf, int64_t subimg_height, int64_t subimg_width,
                                     int64_t subimg_src_origin_x, int64_t subimg_src_origin_y)
{
    MCP3D_ASSERT(subimg_buf)
    mcp3d::TiffDirectoryInfo tiff_info(tif);
    if (tiff_info.is_tiled)
        MCP3D_DOMAIN_ERROR("wrong function call: tiff image is tiled")
    if (subimg_width == 0)
        subimg_width = subimg_height;
    uint32_t img_height = tiff_info.image_height, img_width = tiff_info.image_width, rows_in_strip = tiff_info.rows_in_strip;
    short bits_per_sample = tiff_info.bits_per_sample, samples_per_pixel = tiff_info.samples_per_pixel;
    ValidateTiffSubImageRequest(tiff_info, subimg_height, subimg_width, subimg_src_origin_x, subimg_src_origin_y);
    // tsize_t is return type of TIFFScanlineSize on cluster, instead of tmsize_t
    tsize_t strip_size = TIFFStripSize(const_cast<TIFF*>(tif));
    tdata_t strip_buf = _TIFFmalloc(strip_size);
    if (!strip_buf) MCP3D_BAD_ALLOC("can not allocate memory for tiff strip")
    // scanline does not support random access and require sequentially read
    // through unneeded data. use strips instead
    int64_t subimg_yend = subimg_src_origin_y + subimg_height,
            subimg_xend = subimg_src_origin_x + subimg_width;
    int64_t subimg_sample_offset = 0, strip_sample_offset;
    int64_t copy_y_width = min(subimg_yend, (int64_t)img_height) - subimg_src_origin_y;
    int64_t copy_x_width = min(subimg_xend, (int64_t)img_width) - subimg_src_origin_x;

    // if subimg is partially out of image area, background padding is already
    // done in MImageIO, just make sure to read within bounds
    for (int64_t i = subimg_src_origin_y; i < subimg_src_origin_y + copy_y_width; ++i)
    {
        // offset of the strip within src img
        int64_t strip_offset = i / rows_in_strip;
        // offset of current line offset within strip
        int64_t scanline_offset = i % rows_in_strip;
        // read new strip when previous strip exhausted, or when first entering
        // the loop
        if (scanline_offset == 0 || i == subimg_src_origin_y)
            TIFFReadEncodedStrip(const_cast<TIFF*>(tif), (uint32_t)strip_offset, strip_buf, strip_size);

        strip_sample_offset = samples_per_pixel * (img_width * scanline_offset + subimg_src_origin_x);
        size_t n_bytes = (size_t)samples_per_pixel * copy_x_width * (bits_per_sample / 8);

        if (bits_per_sample == 8)
            memcpy(((uint8_t*)subimg_buf) + subimg_sample_offset, ((uint8_t*)strip_buf) + strip_sample_offset, n_bytes);
        else if (bits_per_sample == 16)
            memcpy(((uint16_t*)subimg_buf) + subimg_sample_offset, ((uint16_t*)strip_buf) + strip_sample_offset, n_bytes);
        else
            memcpy(((int32_t*)subimg_buf) + subimg_sample_offset, ((int32_t*)strip_buf) + strip_sample_offset, n_bytes);

        subimg_sample_offset += samples_per_pixel * subimg_width;
    }
    _TIFFfree(strip_buf);
}

void mcp3d::TransferToSubimg(::TIFF *tif, tdata_t subimg_buf, int64_t subimg_height, int64_t subimg_width,
                             int64_t subimg_src_origin_x, int64_t subimg_src_origin_y)
{
    mcp3d::TiffDirectoryInfo tiff_info(tif);
    if (tiff_info.is_tiled)
        mcp3d::TransferTileToSubimg(tif, subimg_buf, subimg_height, subimg_width, subimg_src_origin_x, subimg_src_origin_y);
    else
        mcp3d::TransferScanlineToSubimg(tif, subimg_buf, subimg_height, subimg_width, subimg_src_origin_x, subimg_src_origin_y);
}

void mcp3d::TransferToSubimg(const string& tif_path, tdata_t subimg_buf, int64_t subimg_height, int64_t subimg_width,
                             int64_t subimg_src_origin_x, int64_t subimg_src_origin_y)
{
    mcp3d::TiffDirectoryInfo tiff_info(tif_path);
    if (tiff_info.is_tiled)
        mcp3d::TransferTileToSubimg(tif_path, subimg_buf, subimg_height, subimg_width, subimg_src_origin_x, subimg_src_origin_y);
    else
        mcp3d::TransferScanlineToSubimg(tif_path, subimg_buf, subimg_height, subimg_width, subimg_src_origin_x, subimg_src_origin_y);
}

