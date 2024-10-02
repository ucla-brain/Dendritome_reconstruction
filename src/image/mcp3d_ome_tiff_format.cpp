//
// Created by muyezhu on 12/2/18.
//
#include <regex>
#include <tiff.h>
#include <tiffio.h>
#include "common/mcp3d_utility.hpp"
#include "image_layout/mcp3d_channel_pyr_slices.hpp"
#include "image_layout/mcp3d_channel_layout.hpp"
#include "image_interface/mcp3d_tiff_utils.hpp"
#include "mcp3d_image_utils.hpp"
#include "image_interface/mcp3d_tiff_io.hpp"
#include "mcp3d_channel_info.hpp"
#include "mcp3d_image.hpp"
#include "mcp3d_image_formats.hpp"
#include "mcp3d_ome_tiff_format.hpp"

using namespace std;


mcp3d::MChannelPyrInfo mcp3d::MOmeTiffFormat::ReadChannelPyrInfo(const mcp3d::MChannelPyrSlices& channel_pyr_slices, int resolution_level)
{
    // assert MChannelPyrSlices has TIFF format
    MCP3D_ASSERT(channel_pyr_slices.file_format() == format_)
    // construct MChannelPyrInfo with channel_pyr_slices
    mcp3d::MChannelPyrInfo channel_pyr_info{channel_pyr_slices};
    // assert resolution_level is consistent with channel_pyr_dir
    MCP3D_ASSERT(mcp3d::MChannelLayout::DirPyrLevel(channel_pyr_info.channel_pyr_dir()) == resolution_level)
    channel_pyr_info.resolution_level_ = resolution_level;
    int n_tiffs = channel_pyr_info.TotalNumberOfImageChunks();
    // because channel_pyr_slices.file_format() was determined to be OMETIFF, n_tiff should be positive
    MCP3D_ASSERT(n_tiffs > 0)
    // MChannelPyrInfo fields based on TiffDirectoryInfo
    mcp3d::TiffDirectoryInfo tiff_info;
    MCP3D_TRY(tiff_info = mcp3d::VerifyTiffSequence(channel_pyr_slices.channel_pyr_dir(), channel_pyr_info.slice_image_names());)
    if (tiff_info.samples_per_pixel > 1)
        MCP3D_MESSAGE("currently only support gray images, rgb iamges will be read as gray")
    if (!mcp3d::GetOmeTiffVolumeHeightWidth(tiff_info.tiff_path, channel_pyr_info.zdim_, channel_pyr_info.ydim_, channel_pyr_info.xdim_))
        MCP3D_RUNTIME_ERROR("can not parse volume zyx dimensions from ome tiff image description")
    channel_pyr_info.chunk_zdim_ = tiff_info.n_directories;
    channel_pyr_info.chunk_ydim_ = tiff_info.image_height;
    channel_pyr_info.chunk_xdim_ = tiff_info.image_width;
    channel_pyr_info.voxel_type_ = mcp3d::TiffBitsPerSampleToVoxelType(tiff_info.bits_per_sample);
    MCP3D_ASSERT(channel_pyr_info.voxel_type_ != mcp3d::VoxelType::UNKNOWN)
    return channel_pyr_info;
}

vector<int> mcp3d::MOmeTiffFormat::OptimalReadXyzDims(const MChannelPyrInfo &pyr_info) const
{
    MCP3D_ASSERT(pyr_info.file_format_ == format_)
    return vector<int>({1, pyr_info.chunk_ydim_, pyr_info.chunk_xdim_});
}

void mcp3d::MOmeTiffFormat::ReadChannelDataImpl(MImage &image, int channel_number)
{
    const MChannelInfo& channel_info = image.image_info().channel_info(channel_number, 0);
    const MImageView& view = image.selected_view();
    int pyr_level = image.selected_view().pyr_level();

    int view_xoffset = view.pyr_level_offsets()[2],
        view_yoffset = view.pyr_level_offsets()[1],
        view_zoffset = view.pyr_level_offsets()[0];
    int view_xextent = view.pyr_level_extents()[2],
        view_yextent = view.pyr_level_extents()[1],
        view_zextent = view.pyr_level_extents()[0];
    int view_xstride = view.pyr_level_strides()[2],
        view_ystride = view.pyr_level_strides()[1],
        view_zstride = view.pyr_level_strides()[0];
    int xchunk_dim = channel_info.xchunk_dim(pyr_level),
        ychunk_dim = channel_info.ychunk_dim(pyr_level),
        zchunk_dim = channel_info.zchunk_dim(pyr_level);

    string tiff_chunk_path;
    // within chunk offsets of source data
    size_t n_bytes_chunk_plane = (size_t)channel_info.xchunk_dim(pyr_level) *
                                 (size_t)channel_info.ychunk_dim(pyr_level) *
                                 (size_t)channel_info.VoxelBytes(pyr_level);
    unique_ptr<uint8_t []> buffer(new (nothrow) uint8_t[n_bytes_chunk_plane]);
    if (!buffer)
        MCP3D_RUNTIME_ERROR("failed to allocate memory for ome tiff plane")
    // traverse yx locations. at each location, read all z planes available in
    // file chunk. this is to reduce file open operations
    int xmax = min(channel_info.xdim(pyr_level), view_xoffset + view_xextent),
        ymax = min(channel_info.ydim(pyr_level), view_yoffset + view_yextent),
        zmax = min(channel_info.zdim(pyr_level), view_zoffset + view_zextent);

    ::TIFF *tiff;
    string ome_tiff_path;
    int x = view_xoffset, y, z;
    int xchunk_offset, ychunk_offset, zchunk_offset, xchunk_id, ychunk_id, zchunk_id;
    int xresidual_stride, yresidual_stride;
    // dst_block tracks read operation's corresponding offsets, extents and
    // strides into view. src_block tracks those of pyr level image
    // the copy operation has offset 0 and extent 1 on z axis in both dst and
    // src blocks due to the ome tif is read one xy plane at a time
    mcp3d::MImageBlock dst_block(vector<int>(3, 0), vector<int>(3, 1), vector<int>(3, 1)),
                       src_block(vector<int>(3, 0), vector<int>(3, 1), vector<int>({1, view_ystride, view_xstride}));
    double time_tiff_open = 0.0, time_tiff_read = 0.0, time_copy = 0.0;
    while (x < xmax)
    {
        xchunk_offset = x % xchunk_dim;
        y = view_yoffset;
        dst_block.offsets_ptr()[1] = 0;
        while (y < ymax)
        {
            ychunk_offset = y % ychunk_dim;
            // retrieve all xy planes at current xy offsets
            z = view_zoffset;
            dst_block.offsets_ptr()[0] = 0;
            dst_block.extents_ptr()[1] = mcp3d::StridedExtent(ychunk_dim - ychunk_offset, view_ystride);
            dst_block.extents_ptr()[2] = mcp3d::StridedExtent(xchunk_dim - xchunk_offset, view_xstride);
            src_block.offsets_ptr()[1] = ychunk_offset;
            src_block.offsets_ptr()[2] = xchunk_offset;
            src_block.extents_ptr()[1] = min(ychunk_dim - ychunk_offset, ymax - y);
            src_block.extents_ptr()[2] = min(xchunk_dim - xchunk_offset, xmax - x);
            while (z < zmax)
            {
                zchunk_id = z / zchunk_dim;
                while (zchunk_id < zmax / zchunk_dim)
                {
                    ome_tiff_path = ImagePath(image, channel_number, pyr_level, z, y, x);
                    INIT_TIMER(0)
                    TIC(0)
                    tiff = TIFFOpen(ome_tiff_path.c_str(), "r");
                    TOC(0)
                    time_tiff_open += ELAPSE(0);
                    // retrieve a file chunk and obtain all data from it
                    // z should be within bounds of image, view and chunk
                    while (z < min(zmax, (zchunk_id + 1) * zchunk_dim))
                    {
                        // read whole z plane
                        zchunk_offset = z % zchunk_dim;
                        INIT_TIMER(1)
                        TIC(1)
                        mcp3d::ReadTiffDirectoryData(tiff, zchunk_offset, buffer);
                        TOC(1)
                        time_tiff_read += ELAPSE(1);
                        // copy to image view
                        INIT_TIMER(2)
                        TIC(2)
                        mcp3d::CopyDataVolume(image.Volume<uint8_t>(channel_number), image.selected_view().xyz_dims(),
                                              buffer.get(), {1, ychunk_dim, xchunk_dim}, image.selected_view().VoxelBytes(), dst_block, src_block);
                        TIC(2)
                        time_copy += ELAPSE(2);
                        // new z level
                        z += view_zstride;
                        // view z offset increments by 1 for next read, while xy offsets remain same
                        dst_block.offsets_ptr()[0] += 1;
                    }
                    // advance zchunk_id
                    TIFFClose(tiff);
                    zchunk_id = z / zchunk_dim;
                }
            }
            // y offset into pyr_level image update
            yresidual_stride = mcp3d::ResidualStride(ychunk_dim - ychunk_offset, view_ystride);
            ychunk_id = y / ychunk_dim;
            if (yresidual_stride <  ychunk_dim)
                y = (ychunk_id + 1) * ychunk_dim + yresidual_stride;
            else
                y += yresidual_stride;
            // y offset into view update
            dst_block.offsets_ptr()[1] += mcp3d::StridedExtent(ychunk_dim - ychunk_offset, view_ystride);
        }
        // x offset into pyr_level image update
        xresidual_stride = mcp3d::ResidualStride(xchunk_dim - xchunk_offset, view_xstride);
        xchunk_id = x / xchunk_dim;
        if (xresidual_stride <  xchunk_dim)
            x = (xchunk_id + 1) * xchunk_dim + xresidual_stride;
        else
            x += xresidual_stride;
        // x offset into view update
        dst_block.offsets_ptr()[2] += mcp3d::StridedExtent(xchunk_dim - xchunk_offset, view_xstride);
    }
    cout << "time to tiff open = " << time_tiff_open << " seconds" << endl;
    cout << "time to tiff read = " << time_tiff_read << " seconds" << endl;
    cout << "time to copy = " << time_copy << " seconds" << endl;
}
