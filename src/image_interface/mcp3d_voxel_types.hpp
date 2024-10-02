//
// Created by muyezhu on 2/11/18.
//
#ifndef MCP3D_MCP3D_SUPPORTED_FORMATS_HPP
#define MCP3D_MCP3D_SUPPORTED_FORMATS_HPP

#include <cstdint>
#include <cfloat>
#include <unordered_set>
#include <map>
#include <hdf5.h>
#include "common/mcp3d_macros.hpp"

namespace mcp3d
{
enum class VoxelType
{
    M8U, M8S, M16U, M16S, M32U, M32S, M64U, M64S, M32F, M64F, UNKNOWN = -1
};

VoxelType StringToVoxelType(const std::string &type_str);

std::string VoxelTypeToString(VoxelType voxel_type);

template <typename VType>
VoxelType TypeToVoxelType()  { return VoxelType::UNKNOWN; }

template<>
inline VoxelType TypeToVoxelType<uint8_t>() { return VoxelType::M8U; }

template<>
inline VoxelType TypeToVoxelType<int8_t>() { return VoxelType::M8S; }

template<>
inline VoxelType TypeToVoxelType<uint16_t>() { return VoxelType::M16U; }

template<>
inline VoxelType TypeToVoxelType<int16_t>() { return VoxelType::M16S; }

template<>
inline VoxelType TypeToVoxelType<uint32_t>() { return VoxelType::M32U; }

template<>
inline VoxelType TypeToVoxelType<int32_t>() { return VoxelType::M32S; }

template<>
inline VoxelType TypeToVoxelType<uint64_t>() { return VoxelType::M64U; }

template<>
inline VoxelType TypeToVoxelType<int64_t>() { return VoxelType::M64S; }

template<>
inline VoxelType TypeToVoxelType<float>() { return VoxelType::M32F; }

template<>
inline VoxelType TypeToVoxelType<double>() { return VoxelType::M64F; }

int VoxelSize(VoxelType voxel_type);

int VoxelTypeToCVType(VoxelType voxel_type, int n_channels = 1);

int MakeCVType(int bytes_per_voxel, int n_channels = 1, bool is_signed = false, bool is_floating_point = false);

template <typename VType>
int TypeToCVType(int n_channels = 1)
{ return mcp3d::VoxelTypeToCVType(mcp3d::TypeToVoxelType<VType>(), n_channels); }

/// invalid argument error is unknown voxel type
hid_t VoxelTypeToH5Type(VoxelType voxel_type);

template <typename VType>
VType VoxelTypeZeroValue();

}


template <typename VType>
VType mcp3d::VoxelTypeZeroValue()
{
    static_assert(std::is_arithmetic<VType>::value, "VType must be arithmetic type");
    if (std::is_same<VType, uint8_t>::value)
        return (uint8_t)0;
    if (std::is_same<VType, int8_t>::value)
        return (int8_t)0;
    if (std::is_same<VType, uint16_t>::value)
        return (uint16_t)0;
    if (std::is_same<VType, int16_t>::value)
        return (int16_t)0;
    if (std::is_same<VType, uint32_t>::value)
        return (uint32_t)0;
    if (std::is_same<VType, int32_t>::value)
        return (int32_t)0;
    if (std::is_same<VType, uint64_t>::value)
        return (uint64_t)0;
    if (std::is_same<VType, int64_t>::value)
        return (int64_t)0;
    if (std::is_same<VType, float>::value)
        return (float)0.0;
    if (std::is_same<VType, double>::value)
        return (double)0.0;
    MCP3D_RUNTIME_ERROR("unknown voxel type")
}

#endif //MCP3D_MCP3D_SUPPORTED_FORMATS_HPP
