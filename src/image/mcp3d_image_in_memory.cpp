//
// Created by muyezhu on 6/13/19.
//

#include "mcp3d_image_in_memory.hpp"

using namespace std;

mcp3d::MImageInMemory::MImageInMemory(const vector<int>& dimensions, VoxelType voxel_type):
        MImageBase(), voxel_type_(voxel_type), is_wrapper_(false), wrapped_data_(std::vector<uint8_t*>{})

{
    MCP3D_ASSERT(voxel_type_ != mcp3d::VoxelType::UNKNOWN)
    set_dimensions(dimensions);
    AllocateSelection();
}

uint8_t* mcp3d::MImageInMemory::data_impl(int volume_index)
{
    if (!is_wrapper_)
    {
        if (volume_index < 0 or volume_index >= (int)data_.size())
            MCP3D_OUT_OF_RANGE("volume index out of range")
        return data_[volume_index].get();
    }
    else
    {
        if (volume_index < 0 or volume_index >= (int)wrapped_data_.size())
            MCP3D_OUT_OF_RANGE("volume index out of range")
        return wrapped_data_[volume_index];
    }
}

const uint8_t* const mcp3d::MImageInMemory::data_impl(int volume_index) const
{
    if (!is_wrapper_)
    {
        if (volume_index < 0 or volume_index >= (int)data_.size())
            MCP3D_OUT_OF_RANGE("volume index out of range")
        return data_[volume_index].get();
    }
    else
    {
        if (volume_index < 0 or volume_index >= (int)wrapped_data_.size())
            MCP3D_OUT_OF_RANGE("volume index out of range")
        return wrapped_data_[volume_index];
    }
}

void mcp3d::MImageInMemory::set_dimensions(const std::vector<int> &dimensions)
{
    MCP3D_ASSERT(dimensions.size() >= (size_t)2 && dimensions.size() <= (size_t)5)
    MCP3D_ASSERT(mcp3d::AllPositive(dimensions))
    dimensions_ = mcp3d::To5D(dimensions);
    // supports a single time point
    MCP3D_ASSERT(dimensions_[0] == 1)
}

void mcp3d::MImageInMemory::AllocateSelection()
{
    MCP3D_ASSERT(!is_wrapper_)
    MCP3D_ASSERT(!empty())
    data_.clear();
    auto n_volume_bytes = mcp3d::ReduceProdSeq<size_t>(xyz_dims());
    for (int i = 0; i < dimensions_[1]; ++i)
    {
        data_.push_back(make_unique<uint8_t[]>((size_t)mcp3d::VoxelSize(voxel_type_) * n_volume_bytes));
        MCP3D_ASSERT(data_[i])
    }
}

void mcp3d::MImageInMemory::Clear()
{
    dimensions_.clear();
    voxel_type_ = mcp3d::VoxelType::UNKNOWN;
    is_wrapper_ = false;
    data_.clear();
    wrapped_data_.clear();
}
