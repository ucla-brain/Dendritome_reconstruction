//
// Created by muyezhu on 3/11/19.
//
#include <algorithm>
#include "common/mcp3d_common.hpp"
#include "mcp3d_hdf5_utils.hpp"
#include "mcp3d_imaris_util.hpp"

using namespace std;

vector<string> mcp3d::ImarisResolutionNames(const string &imaris_path)
{
    hid_t imaris_id = mcp3d::Hdf5Handle(imaris_path);
    hid_t dataset_group_id = H5Gopen(imaris_id, "DataSet", H5P_DEFAULT);
    MCP3D_ASSERT(dataset_group_id >= 0)
    vector<string> object_names(mcp3d::ObjectNamesInGroup(dataset_group_id));
    H5Gclose(dataset_group_id);
    H5Fclose(imaris_id);
    return object_names;
}

vector<string> mcp3d::ImarisResolutionChannelNames(hid_t imaris_id, int resolution_level, int time)
{
    hid_t time_point_id = mcp3d::ImarisTimePointHandle(imaris_id, resolution_level, time);
    vector<string> object_names(mcp3d::ObjectNamesInGroup(time_point_id));
    sort(object_names.begin(), object_names.end());
    H5Gclose(time_point_id);
    return object_names;
}

vector<string> mcp3d::ImarisResolutionChannelNames(const string &imaris_path, int resolution_level, int time)
{
    hid_t imaris_id = mcp3d::Hdf5Handle(imaris_path);
    vector<string> channel_names(mcp3d::ImarisResolutionChannelNames(imaris_id, resolution_level, time));
    H5Fclose(imaris_id);
    return channel_names;
}

vector<int> mcp3d::ImarisChannelImageXyzSizes(hid_t imaris_id, int resolution_level, int channel_number, int time)
{
    hid_t channel_id = mcp3d::ImarisChannelHandle(imaris_id, resolution_level, channel_number, time);
    int zdim, ydim, xdim;
    zdim = stoi(mcp3d::Hdf5AttributeString(channel_id, "ImageSizeZ"));
    ydim = stoi(mcp3d::Hdf5AttributeString(channel_id, "ImageSizeY"));
    xdim = stoi(mcp3d::Hdf5AttributeString(channel_id, "ImageSizeX"));
    H5Gclose(channel_id);
    return vector<int>({zdim, ydim, xdim});
}

vector<int> mcp3d::ImarisResolutionImageXyzSizes(hid_t imaris_id, int resolution_level, int time)
{
    vector<string> channel_names = mcp3d::ImarisResolutionChannelNames(imaris_id, resolution_level, time);
    vector<int> xyz_sizes = mcp3d::ImarisChannelImageXyzSizes(imaris_id, resolution_level, 0, time);
    for(int i = 1; i < (int)channel_names.size(); ++i)
        if (xyz_sizes != mcp3d::ImarisChannelImageXyzSizes(imaris_id, resolution_level, i, time))
            MCP3D_RUNTIME_ERROR("imaris files have inconsistent image XYZ sizes across channels at resolution level " + to_string(resolution_level))
    return xyz_sizes;
}

vector<int> mcp3d::ImarisChannelChunkXyzDims(hid_t imaris_id, int resolution_level, int channel_number, int time)
{
    hid_t dataset_id = mcp3d::ImarisDatasetHandle(imaris_id, resolution_level, channel_number, time);
    vector<int> chunk_dims = mcp3d::Hdf5DatasetChunkDimensions(dataset_id);
    H5Dclose(dataset_id);
    return chunk_dims;
}

vector<int> mcp3d::ImarisResolutionChunkXyzDims(hid_t imaris_id, int resolution_level, int time)
{
    vector<string> channel_names = mcp3d::ImarisResolutionChannelNames(imaris_id, resolution_level, time);
    vector<int> chunk_dims = mcp3d::ImarisChannelChunkXyzDims(imaris_id, resolution_level, 0, time);
    for(int i = 1; i < (int)channel_names.size(); ++i)
        if (chunk_dims != mcp3d::ImarisChannelChunkXyzDims(imaris_id, resolution_level, i, time))
        MCP3D_RUNTIME_ERROR("imaris files have inconsistent chunk dimensions across channels at resolution level " + to_string(resolution_level))
    return chunk_dims;
}

mcp3d::VoxelType mcp3d::ImarisChannelVoxelType(hid_t imaris_id, int resolution_level, int channel_number, int time)
{
    hid_t dataset_id = mcp3d::ImarisDatasetHandle(imaris_id, resolution_level, channel_number, time);
    return mcp3d::Hdf5DatasetVoxelType(dataset_id);
}

mcp3d::VoxelType mcp3d::ImarisResolutionVoxelType(hid_t imaris_id, int resolution_level, int time)
{
    vector<string> channel_names = mcp3d::ImarisResolutionChannelNames(imaris_id, resolution_level, time);
    mcp3d::VoxelType voxel_type = mcp3d::ImarisChannelVoxelType(imaris_id, resolution_level, 0, time);
    MCP3D_ASSERT(voxel_type != mcp3d::VoxelType::UNKNOWN)
    for(int i = 1; i < (int)channel_names.size(); ++i)
        if (voxel_type != mcp3d::ImarisChannelVoxelType(imaris_id, resolution_level, i, time))
            MCP3D_RUNTIME_ERROR("imaris files have inconsistent voxel datatype across channels at resolution level " + to_string(resolution_level))
    return voxel_type;
}

hid_t mcp3d::ImarisResolutionHandle(hid_t imaris_id, int resolution_level)
{
    MCP3D_ASSERT(imaris_id >= 0)
    MCP3D_ASSERT(resolution_level >= 0)
    string res_group_path("/DataSet/ResolutionLevel ");
    res_group_path.append(to_string(resolution_level));
    hid_t resolution_id = H5Gopen(imaris_id, res_group_path.c_str(), H5P_DEFAULT);
    if (resolution_id < 0)
        MCP3D_RUNTIME_ERROR("can not open group at resolution level " + to_string(resolution_level))
    return resolution_id;
}

hid_t mcp3d::ImarisTimePointHandle(hid_t imaris_id, int resolution_level, int time)
{
    MCP3D_ASSERT(time >= 0)
    hid_t resolution_id = mcp3d::ImarisResolutionHandle(imaris_id, resolution_level);
    string time_point_path("TimePoint ");
    time_point_path.append(to_string(time));
    hid_t time_point_id = H5Gopen(resolution_id, time_point_path.c_str(), H5P_DEFAULT);
    if (time_point_id < 0)
        MCP3D_RUNTIME_ERROR("can not open group at resolution level " + to_string(resolution_level) + " time point " + to_string(time))
    H5Gclose(resolution_id);
    return time_point_id;
}

hid_t mcp3d::ImarisChannelHandle(hid_t imaris_id, int resolution_level, int channel_number, int time)
{
    MCP3D_ASSERT(channel_number >= 0)
    hid_t time_point_id = mcp3d::ImarisTimePointHandle(imaris_id, resolution_level, time);
    string channel_name("Channel ");
    channel_name.append(to_string(channel_number));
    hid_t channel_id = H5Gopen(time_point_id, channel_name.c_str(), H5P_DEFAULT);
    if (channel_id < 0)
        MCP3D_RUNTIME_ERROR("can not open group at resolution level " + to_string(resolution_level) + " time point " +
                            to_string(time) + " channel number " + to_string(channel_number))
    H5Gclose(time_point_id);
    return channel_id;
}

hid_t mcp3d::ImarisDatasetHandle(hid_t imaris_id, int resolution_level, int channel_number, int time)
{
    hid_t channel_id = mcp3d::ImarisChannelHandle(imaris_id, resolution_level, channel_number, time);
    hid_t data_id = H5Dopen(channel_id, "Data", H5P_DEFAULT);
    if (data_id < 0)
        MCP3D_RUNTIME_ERROR("can not open dataset at resolution level " + to_string(resolution_level) + " time point " +
                            to_string(time) + " channel number " + to_string(channel_number))
    H5Gclose(channel_id);
    return data_id;
}



