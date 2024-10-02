//
// Created by muyezhu on 4/5/19.
//

#ifndef MCP3D_MCP3D_IMAGE_BASE_HPP
#define MCP3D_MCP3D_IMAGE_BASE_HPP

#include <cstdint>
#include <memory>
#include "image_interface/mcp3d_voxel_types.hpp"

namespace mcp3d
{

class MImageBase
{
public:
    MImageBase(): data_(std::vector<std::unique_ptr<uint8_t []>>{}),
                  older_data_(std::vector<std::vector<std::unique_ptr<uint8_t []>>>{}) {}

    virtual ~MImageBase()= default;

    virtual bool empty() const = 0;

    virtual void ClearLoaded() = 0;

    /// cast data in data_ vector
    void Cast(VoxelType target_type);

    void Normalize(VoxelType target_type = VoxelType::UNKNOWN);

    /// return false if neither has loaded data
    /// return true if images have same loaded_dims(), and each volume holds identical data
    bool HasEqualData(const MImageBase &other);

    virtual std::vector<int> loaded_dims() const = 0;

    virtual std::vector<int> loaded_xyz_dims() const = 0;

    virtual int loaded_n_channels() const = 0;

    virtual VoxelType loaded_voxel_type() const = 0;

    // instance getters
    /// returns a uint8_t pointer to data_[volume_index] (or wrapped_data_ if instance is wrapper). note the use of volume_index
    /// this function is intended to make the base class able to perform common array based operations such as cast
    uint8_t* data(int volume_index = 0)
    { return data_impl(volume_index); }

    const uint8_t* const data(int volume_index = 0) const
    { return data_impl(volume_index); }

    friend class MImageIO;

protected:
    int64_t VoxelAddress(int c, int z, int y, int x, int t = 0) const;

    virtual uint8_t* data_impl(int volume_index) = 0;

    virtual const uint8_t* const data_impl(int volume_index) const = 0;

    /// each std::unique_ptr<uint8_t []> manages a single zyx volue
    /// data_.size() is equal to number of loaded channels
    std::vector<std::unique_ptr<uint8_t []>> data_;
    std::vector<std::vector<std::unique_ptr<uint8_t []>>> older_data_;

private:
    virtual void AllocateSelection() = 0;
};

}


#endif //MCP3D_MCP3D_IMAGE_BASE_HPP
