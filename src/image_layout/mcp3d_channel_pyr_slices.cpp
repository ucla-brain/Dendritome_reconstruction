//
// Created by muyezhu on 3/17/19.
//
#include <algorithm>
#include "image_interface/mcp3d_imaris_util.hpp"
#include "mcp3d_channel_layout.hpp"
#include "mcp3d_channel_pyr_slices.hpp"

using namespace std;

bool mcp3d::MChannelPyrSlices::operator==(const mcp3d::MChannelPyrSlices& other) const
{
    return channel_pyr_dir_ == other.channel_pyr_dir_ &&
           slice_dims_ == other.slice_dims_ &&
           slice_numbers_ == other.slice_numbers_ &&
           file_format_ == other.file_format_;
}

void mcp3d::MChannelPyrSlices::Refresh()
{
    if (!mcp3d::IsDir(channel_pyr_dir_))
    {
        Clear();
        return;
    }
    ReadChannelPyrSlices();
    ReadChannelPyrFormat();
}

void mcp3d::MChannelPyrSlices::Clear()
{
    channel_pyr_dir_.clear();
    slice_dims_.clear();
    slice_numbers_.clear();
    file_format_ = mcp3d::FileFormat::UNKNOWN;
}

string mcp3d::MChannelPyrSlices::AxisSliceNameFromSliceID(mcp3d::ChannelAxis axis, int slice_id) const
{
    if (empty())
        return "";
    MCP3D_ASSERT(slice_id >= 0 && slice_id < slice_numbers_.at(axis))
    int slice_dim = slice_dims_.at(axis);
    return mcp3d::MChannelPyrSlices::AxisSliceName(axis, slice_id, slice_dim);
}

string mcp3d::MChannelPyrSlices::AxisSliceNameFromCoordinate(ChannelAxis axis, int coordinate) const
{
    if (empty())
        return "";
    MCP3D_ASSERT(coordinate >= 0)
    if (slice_dims_.at(axis) != mcp3d::MAX_AXIS_SLICE_DIM)
        MCP3D_ASSERT(coordinate < slice_dims_.at(axis) * slice_numbers_.at(axis))
    int slice_id = slice_dims_.at(axis) == mcp3d::MAX_AXIS_SLICE_DIM ? 0 : coordinate / slice_dims_.at(axis);
    return mcp3d::MChannelPyrSlices::AxisSliceName(axis, slice_id, slice_dims_.at(axis));
}

string mcp3d::MChannelPyrSlices::AxisSliceName(mcp3d::ChannelAxis axis, int slice_id, int slice_dim)
{
    if (slice_dim == mcp3d::MAX_AXIS_SLICE_DIM)
        return "";
    MCP3D_ASSERT(slice_dim >= 0 && (int)to_string(slice_dim).size() <= mcp3d::SLICE_DIM_WIDTH)
    MCP3D_ASSERT(slice_id >= 0 && (int)to_string(slice_id).size() <= mcp3d::SLICE_ID_WIDTH)
    string slice_name(ChannelAxisStr(axis));
    slice_name.append(mcp3d::PadNumStr(slice_id, mcp3d::SLICE_ID_WIDTH));
    slice_name.append("_");
    string size_string = mcp3d::PadNumStr(slice_dim, mcp3d::SLICE_DIM_WIDTH);
    slice_name.append(size_string);
    return slice_name;
}

string mcp3d::MChannelPyrSlices::SliceNameFromSliceIDs(int zslice_id, int yslice_id, int xslice_id) const
{
    return mcp3d::JoinPath(AxisSliceNameFromSliceID(mcp3d::ChannelAxis::Z, zslice_id),
                           AxisSliceNameFromSliceID(mcp3d::ChannelAxis::Y, yslice_id),
                           AxisSliceNameFromSliceID(mcp3d::ChannelAxis::X, xslice_id));
}

string mcp3d::MChannelPyrSlices::SliceNameFromCoordinates(int z_coordinate, int y_coordinate, int x_coordinate) const
{
    return mcp3d::JoinPath(AxisSliceNameFromCoordinate(mcp3d::ChannelAxis::Z, z_coordinate),
                           AxisSliceNameFromCoordinate(mcp3d::ChannelAxis::Y, y_coordinate),
                           AxisSliceNameFromCoordinate(mcp3d::ChannelAxis::X, x_coordinate));
}

vector<string> mcp3d::MChannelPyrSlices::SliceNames() const
{
    vector<string> names;
    if (empty())
        return names;
    for (int i = 0; i < slice_numbers_.at(mcp3d::ChannelAxis::Z); ++i)
        for (int j = 0; j < slice_numbers_.at(mcp3d::ChannelAxis::Y); ++j)
            for (int k = 0; k < slice_numbers_.at(mcp3d::ChannelAxis::X); ++k)
            {
                string name = mcp3d::JoinPath(AxisSliceNameFromSliceID(mcp3d::ChannelAxis::Z, i),
                                              AxisSliceNameFromSliceID(mcp3d::ChannelAxis::Y, j),
                                              AxisSliceNameFromSliceID(mcp3d::ChannelAxis::X, k));
                names.push_back(name);
            }
    return names;
}

void mcp3d::MChannelPyrSlices::ReadChannelPyrSlices()
{
    if (empty())
        return;
    slice_dims_.clear();
    slice_numbers_.clear();
    file_format_ = mcp3d::FileFormat::UNKNOWN;
    vector<string> zslice_names, yslice_names, xslice_names;
    unordered_map<string, unordered_set<int>> slice_ids {};
    ReadChannelPyrSlicesTree(channel_pyr_dir_, slice_ids);
    // update slice_numbers_ and if needed slice_dims_
    mcp3d::ChannelAxis axis;
    for (const auto& axis_str: vector<string>({"z", "y", "x"}))
    {
        axis = mcp3d::StrToChannelAxis(axis_str);
        if (slice_ids.find(axis_str) != slice_ids.end())
            slice_numbers_[axis] = (int)slice_ids[axis_str].size();
        else
        {
            slice_numbers_[axis] = 1;
            slice_dims_[axis] = mcp3d::MAX_AXIS_SLICE_DIM;
        }
        if ((slice_numbers_[axis] == 1) ^ (slice_dims_[axis] == mcp3d::MAX_AXIS_SLICE_DIM))
            MCP3D_RUNTIME_ERROR("axis with a single slice should always have mcp3d::MAX_AXIS_SLICE_DIM as slice dimension, and vice versa")
    }
    for (const auto& slice_name: SliceNames())
    {
        string slice_dir{mcp3d::JoinPath(channel_pyr_dir_, slice_name)};
        if (!mcp3d::IsDir(slice_dir))
            MCP3D_RUNTIME_ERROR("expected slice directory " + slice_dir + " not found")
    }
}

void mcp3d::MChannelPyrSlices::ReadChannelPyrSlicesTree(const string &input_dir, unordered_map<std::string, unordered_set<int>>& slice_ids)
{
    // process input_dir if its a folder slice dir
    string input_axis_str, basename = mcp3d::Basename(input_dir);
    smatch sm;
    if (regex_match(basename, sm, mcp3d::MChannelPyrSlices::AxisSliceNamePattern()))
    {
        input_axis_str = sm.str(1);
        mcp3d::ChannelAxis input_axis = mcp3d::StrToChannelAxis(input_axis_str);
        // place slide_id in slice_ids[axis]
        int slice_id = stoi(sm.str(2));
        slice_ids[input_axis_str].insert(slice_id);
        // assert that parsed slice_dim for each axis are consistent
        int slice_dim = stoi(sm.str(3));
        if (slice_dims_.find(input_axis) == slice_dims_.end())
            slice_dims_[input_axis] = slice_dim;
        else
            MCP3D_ASSERT(slice_dims_[input_axis] == slice_dim)
    }

    // if slice folders found under input_dir, recursively process them
    vector<string> slice_dirs;
    string common_axis_str, dir_axis_str;
    for (const auto& slice_name: mcp3d::DirsInDir(input_dir, false, false))
    {
        if (regex_match(slice_name, sm, mcp3d::MChannelPyrSlices::AxisSliceNamePattern()))
        {
            dir_axis_str = sm.str(1);
            // check that folder slices under input_dir belong to the same axis
            if (common_axis_str.empty())
                common_axis_str = dir_axis_str;
            else if (common_axis_str != dir_axis_str)
                MCP3D_RUNTIME_ERROR("folder slices under " + input_dir + " do not belong to the same axis")
            slice_dirs.push_back(mcp3d::JoinPath(input_dir, slice_name));
        }
    }
    // check that the common axis follows zyx ordering compared to input_dir's axis
    if (!input_axis_str.empty() && !common_axis_str.empty())
        if (input_axis_str <= common_axis_str)
            MCP3D_RUNTIME_ERROR("folder slices under " + input_dir + " do not follow the expected zyx axes ordering")

    // recursive call
    for (const auto& slice_dir: slice_dirs)
        ReadChannelPyrSlicesTree(slice_dir, slice_ids);
}

void mcp3d::MChannelPyrSlices::ReadChannelPyrFormat()
{
    if (empty())
        return;
    for (const auto& format: mcp3d::FileFormatSearchOrder())
    {
        if (!ChannelPyrDirCompatibleWithFileFormat(format))
            continue;
        for (const auto& slice_name: SliceNames())
            if (mcp3d::DirContainsFileFormat(mcp3d::JoinPath(channel_pyr_dir_, slice_name), format))
            {
                file_format_ = format;
                return;
            }
    }
    file_format_ = mcp3d::FileFormat::UNKNOWN;
}

bool mcp3d::MChannelPyrSlices::ChannelPyrDirCompatibleWithFileFormat(FileFormat file_format)
{
    if (file_format == mcp3d::FileFormat::TIFF)
        return ChannelPyrDirCompatibleWithTiffFormat();
    else if (file_format == mcp3d::FileFormat::OMETIFF)
        return ChannelPyrDirCompatibleWithOmeTiffFormat();
    else if (file_format == mcp3d::FileFormat::IMARIS)
        return ChannelPyrDirCompatibleWithImarisFormat();
    else if (file_format == mcp3d::FileFormat::HDF5)
        return ChannelPyrDirCompatibleWithHdf5Format();
    else
        return false;
}

void mcp3d::MChannelPyrSlices::set_channel_pyr_dir(const string &channel_pyr_dir)
{
    if (!mcp3d::IsDir(channel_pyr_dir))
        MCP3D_OS_ERROR(channel_pyr_dir + " is not a directory")
    channel_pyr_dir_ = channel_pyr_dir;
    Refresh();
}
