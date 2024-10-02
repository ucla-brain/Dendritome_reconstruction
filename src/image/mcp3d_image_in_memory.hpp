//
// Created by muyezhu on 6/13/19.
//

#ifndef MCP3D_MCP3D_IMAGE_IN_MEMORY_HPP
#define MCP3D_MCP3D_IMAGE_IN_MEMORY_HPP

#include "common/mcp3d_common.hpp"
#include "mcp3d_image_base.hpp"

namespace mcp3d
{

class MImageInMemory: public MImageBase
{
public:
    /// default constructed instances are empty (empty() returns true)
    MImageInMemory(): dimensions_(5, 0), voxel_type_(VoxelType::UNKNOWN), is_wrapper_(false), wrapped_data_(std::vector<uint8_t*>{}) {}

    /// MImageInMemory instance has a single evel 0
    /// dimensions must have at least 2 and at most 5 positive integers, and will be expanded
    /// to 5d vector representing tczyx dimenions
    /// allocate memory according to dimensions
    explicit MImageInMemory(const std::vector<int>& dimensions, VoxelType voxel_type = VoxelType::M16U);

    MImageInMemory(const MImageInMemory& other) = delete;

    bool empty() const override
    { return AllPositive(dimensions_); }

    template <typename VType = uint8_t>
    VType* Volume(int channel_number = 0, int time = 0);

    /// read only exposure to data pointer
    template <typename VType = uint8_t>
    const VType* ConstVolume(int channel_number = 0, int time = 0) const;

    template <typename VType = uint8_t>
    VType* Plane(int channel_number, int view_z, int time = 0);

    /// the coordinates here are local to voxels in loaded image view.
    /// first element in first array of the data_ vector is indexed (0, 0, 0, 0, 0)
    template <typename VType = uint8_t>
    void SetVoxel(int c, int z, int y, int x, VType val, int t = 0)
    { Volume<VType>(c, t)[VoxelAddress(c, z, y, x)] = val; }

    template <typename VType = uint8_t>
    VType& operator() (int c, int z, int y, int x, int t = 0)
    { return Volume<VType>(c, t)[VoxelAddress(c, z, y, x)]; }

    template <typename VType = uint8_t>
    const VType& At(int c, int z, int y, int x, int t = 0) const
    { return ConstVolume<VType>(c, t)[VoxelAddress(c, z, y, x)]; }

    // clears data, wrapped_data, dimensions_, set voxel_type_ to unknwon and
    // is_wrapper_ to false
    void Clear();

    bool is_wrapper() const
    { return is_wrapper_; }

    /// does not take ownership of data. borrows raw pointer instead
    template <typename VType = uint8_t>
    void WrapData(std::vector<VType*> src_ptrs, const std::vector<int>& src_dims);

    int zdim() const
    { return dimensions_[2]; }

    int ydim() const
    { return dimensions_[3]; }

    int xdim() const
    { return dimensions_[4]; }

    std::vector<int> xyz_dims() const
    { return std::vector<int>({dimensions_[2], dimensions_[3], dimensions_[4]}); }

    int n_channels() const
    { return dimensions_[1]; }

    int n_times() const
    { return dimensions_[0]; }

    int n_volumes() const
    { return n_channels() * n_times(); }

    std::vector<int> loaded_dims() const override
    { return dimensions_; }

    std::vector<int> loaded_xyz_dims() const override
    { return xyz_dims(); }

    int loaded_n_channels() const override
    { return dimensions_[1]; }

    VoxelType loaded_voxel_type() const override
    { return empty() ? VoxelType::UNKNOWN : voxel_type_; }

private:
    /// asserts !is_wrapper_. wrappers operates on pre-allocated data
    void AllocateSelection() override;

    uint8_t* data_impl(int volume_index) override;

    const uint8_t* const data_impl(int volume_index) const override;

    void set_dimensions(const std::vector<int>& dimensions);

    std::vector<int> dimensions_;
    VoxelType voxel_type_;
    bool is_wrapper_;
    std::vector<uint8_t*> wrapped_data_;

};

}

template <typename VType>
VType* mcp3d::MImageInMemory::Volume(int channel_number, int time)
{
    MCP3D_ASSERT(!empty())
    MCP3D_ASSERT(channel_number >= 0 && channel_number < n_channels())
    MCP3D_ASSERT(time >= 0 && time < n_times())
    if (!is_wrapper_)
        return reinterpret_cast<VType*>(data_[channel_number].get());
    else
        return reinterpret_cast<VType*>(wrapped_data_[channel_number]);
}

template <typename VType>
const VType* mcp3d::MImageInMemory::ConstVolume(int channel_number, int time) const
{
    MCP3D_ASSERT(!empty())
    MCP3D_ASSERT(channel_number >= 0 && channel_number < n_channels())
    MCP3D_ASSERT(time >= 0 && time < n_times())
    if (!is_wrapper_)
        return reinterpret_cast<const VType*>(data_[channel_number].get());
    else
        return reinterpret_cast<const VType*>(wrapped_data_[channel_number]);
}


template <typename VType>
VType* mcp3d::MImageInMemory::Plane(int channel_number, int view_z, int time)
{
    long n_plane_voxels = (long)dimensions_[3] * (long)dimensions_[4];
    return Volume<VType>(channel_number, time) + n_plane_voxels * (long)view_z;
}

template <typename VType>
void mcp3d::MImageInMemory::WrapData(std::vector<VType*> src_ptrs, const std::vector<int>& src_dims)
{
    Clear();
    voxel_type_ = mcp3d::TypeToVoxelType<VType>();
    is_wrapper_ = true;
    set_dimensions(src_dims);
    for (const auto& src_ptr: src_ptrs)
        wrapped_data_.push_back(reinterpret_cast<uint8_t*>(src_ptr));
}

#endif //MCP3D_MCP3D_IMAGE_IN_MEMORY_HPP
