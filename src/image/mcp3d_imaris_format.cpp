// Created by muyezhu on 1/30/19.
#include <cstring>
#include "image_interface/mcp3d_voxel_types.hpp"
#include "image_interface/mcp3d_hdf5_utils.hpp"
#include "image_interface/mcp3d_imaris_util.hpp"
#include "image_layout/mcp3d_channel_pyr_slices.hpp"
#include "image_layout/mcp3d_channel_layout.hpp"
#include "mcp3d_channel_info.hpp"
#include "mcp3d_image.hpp"
#include "mcp3d_imaris_format.hpp"

using namespace std;


mcp3d::ImarisResolutionLevelInfo::ImarisResolutionLevelInfo(const string &imaris_path, int resolution_level, int time):
        imaris_path_(imaris_path), time_(time), resolution_level_(resolution_level)
{
    hid_t imaris_id = mcp3d::Hdf5Handle(imaris_path);
    image_xyz_sizes_ = mcp3d::ImarisResolutionImageXyzSizes(imaris_id, resolution_level_, time_);
    chunk_xyz_dims_ = mcp3d::ImarisResolutionChunkXyzDims(imaris_id, resolution_level_, time_);
    voxel_type_ = mcp3d::ImarisResolutionVoxelType(imaris_id, resolution_level_, time_);
    H5Fclose(imaris_id);
}

mcp3d::MChannelPyrInfo mcp3d::MImarisFormat::ReadChannelPyrInfo(const mcp3d::MChannelPyrSlices &channel_pyr_slices, int resolution_level)
{
    // assert MChannelPyrSlices has IMARIS format
    MCP3D_ASSERT(channel_pyr_slices.file_format() == format_)
    // construct MChannelPyrInfo with channel_pyr_slices
    mcp3d::MChannelPyrInfo channel_pyr_info{channel_pyr_slices};
    MCP3D_ASSERT(mcp3d::MChannelLayout::DirPyrLevel(channel_pyr_info.channel_pyr_dir()) == 0)
    string imaris_path = mcp3d::StitchedImarisPathInDir(channel_pyr_info.channel_pyr_dir());
    // assert resolution level exist in imaris file
    MCP3D_ASSERT(resolution_level >= 0 && resolution_level < mcp3d::NumberOfImarisResolutions(imaris_path))
    mcp3d::ImarisResolutionLevelInfo imaris_info(imaris_path, resolution_level, 0);
    channel_pyr_info.resolution_level_ = resolution_level;
    channel_pyr_info.zdim_ = imaris_info.image_xyz_sizes()[0];
    channel_pyr_info.ydim_ = imaris_info.image_xyz_sizes()[1];
    channel_pyr_info.xdim_ = imaris_info.image_xyz_sizes()[2];
    // single imaris physical file
    channel_pyr_info.chunk_zdim_ = channel_pyr_info.zdim_;
    channel_pyr_info.chunk_ydim_ = channel_pyr_info.ydim_;
    channel_pyr_info.chunk_xdim_ = channel_pyr_info.xdim_;
    channel_pyr_info.voxel_type_ = imaris_info.voxel_type();
    return channel_pyr_info;
}

vector<int> mcp3d::MImarisFormat::OptimalReadXyzDims(const mcp3d::MChannelPyrInfo &pyr_info) const
{
    MCP3D_ASSERT(pyr_info.file_format_ == format_)
    string imaris_path(pyr_info.slice_image_names_.at("")[0]);
    mcp3d::ImarisResolutionLevelInfo imaris_info(imaris_path, pyr_info.resolution_level_, 0);
    if (imaris_info.is_chunked())
        return imaris_info.chunk_xyz_dims();
    else
        return imaris_info.image_xyz_sizes();
}

void mcp3d::MImarisFormat::ReadChannelDataImpl(MImage &image, int channel_number)
{
    const mcp3d::MImageView& view = image.selected_view();
    int pyr_level = view.pyr_level();
    string image_path(ImagePath(image, channel_number, pyr_level, 0, 0, 0));
    hid_t hdf5_id = mcp3d::Hdf5Handle(image_path);
    // assuming single time point
    hid_t dataset_id = mcp3d::ImarisDatasetHandle(hdf5_id, pyr_level, channel_number, 0);
    hsize_t file_ds_starts[3] = {(hsize_t)view.pyr_level_offsets()[0], (hsize_t)view.pyr_level_offsets()[1], (hsize_t)view.pyr_level_offsets()[2]};
    hsize_t file_ds_strides[3] = {(hsize_t)view.pyr_level_strides()[0], (hsize_t)view.pyr_level_strides()[1], (hsize_t)view.pyr_level_strides()[2]};
    vector<int> counts = view.xyz_dims_in_bound();
    hsize_t ds_counts[3] = {(hsize_t)counts[0], (hsize_t)counts[1], (hsize_t)counts[2]};
    hsize_t mem_ds_dims[3] = {(hsize_t)view.zdim(), (hsize_t)view.ydim(), (hsize_t)view.xdim()};
    hsize_t mem_ds_starts[3] = {0, 0, 0};
    hid_t mem_type_id = mcp3d::VoxelTypeToH5Type(view.voxel_type());
    mcp3d::ReadHdf5Dataset(dataset_id, file_ds_starts, file_ds_strides, ds_counts, mem_ds_dims, mem_ds_starts,
                           mem_type_id, image.Volume<uint8_t>(channel_number));
    H5Dclose(dataset_id);
    H5Fclose(hdf5_id);
}