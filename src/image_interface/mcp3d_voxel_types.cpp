//
// Created by muyezhu on 2/27/18.
//
#include <opencv2/core/core.hpp>
#include "common/mcp3d_common.hpp"
#include "mcp3d_voxel_types.hpp"

using namespace std;

mcp3d::VoxelType mcp3d::StringToVoxelType(const string &type_str)
{
    if (type_str == "M8U")
        return VoxelType::M8U;
    else if (type_str == "M8S")
        return VoxelType::M8S;
    else if (type_str == "M16U")
        return VoxelType::M16U;
    else if (type_str == "M16S")
        return VoxelType::M16S;
    else if (type_str == "M32U")
        return VoxelType::M32U;
    else if (type_str == "M32S")
        return VoxelType::M32S;
    else if (type_str == "M64U")
        return VoxelType::M64U;
    else if (type_str == "M64S")
        return VoxelType::M64S;
    else if (type_str == "M32F")
        return VoxelType::M32F;
    else if (type_str == "M64F")
        return VoxelType::M64F;
    else
        return VoxelType::UNKNOWN;
}

string mcp3d::VoxelTypeToString(VoxelType voxel_type)
{
    if (voxel_type == VoxelType::M8U)
        return "M8U";
    else if (voxel_type == VoxelType::M8S)
        return "M8S";
    else if (voxel_type == VoxelType::M16U)
        return "M16U";
    else if (voxel_type == VoxelType::M16S)
        return "M16S";
    else if (voxel_type == VoxelType::M32U)
        return "M32U";
    else if (voxel_type == VoxelType::M32S)
        return "M32S";
    else if (voxel_type == VoxelType::M64U)
        return "M64U";
    else if (voxel_type == VoxelType::M64S)
        return "M64S";
    else if (voxel_type == VoxelType::M32F)
        return "M32F";
    else if (voxel_type == VoxelType::M64F)
        return "M64F";
    else
        return "UNKNOWN";
}

int mcp3d::VoxelSize(VoxelType voxel_type)
{
    if (voxel_type == VoxelType::M8U || voxel_type == VoxelType::M8S)
        return 1;
    else if (voxel_type == VoxelType::M16U || voxel_type == VoxelType::M16S)
        return 2;
    else if (voxel_type == VoxelType::M32U || voxel_type == VoxelType::M32S || voxel_type == VoxelType::M32F)
        return 4;
    else if (voxel_type == VoxelType::M64U || voxel_type == VoxelType::M64S || voxel_type == VoxelType::M64F)
        return 8;
    else
        return 0;
}

int mcp3d::VoxelTypeToCVType(VoxelType voxel_type, int n_channels)
{
    MCP3D_ASSERT(n_channels >= 1)
    switch (voxel_type)
    {
        case VoxelType::M8U:
            return CV_8UC(n_channels);
        case VoxelType::M8S:
            return CV_8SC(n_channels);
        case VoxelType::M16U:
            return CV_16UC(n_channels);
        case VoxelType::M16S:
            return CV_16SC(n_channels);
        case VoxelType::M32U:
            cout << "warning: opencv does not have CV_32U, returning CV_32S" << endl;
            return CV_32SC(n_channels);
        case VoxelType::M32S:
            return CV_32SC(n_channels);
        case VoxelType::M64U:
            MCP3D_DOMAIN_ERROR("opencv does not have 64bit integer CvType")
        case VoxelType::M64S:
            MCP3D_DOMAIN_ERROR("opencv does not have 64bit integer CvType")
        case VoxelType::M32F:
            return CV_32FC(n_channels);
        case VoxelType::M64F:
            return CV_64FC(n_channels);
        default:
            MCP3D_INVALID_ARGUMENT("unknown VoxelType value")
    }
}

int mcp3d::MakeCVType(int bytes_per_voxel, int n_channels, bool is_signed, bool is_floating_point)
{
    MCP3D_ASSERT(bytes_per_voxel == 1 || bytes_per_voxel == 2 ||
                 bytes_per_voxel == 4 || bytes_per_voxel == 8)
    MCP3D_ASSERT(n_channels >= 1)
    switch (bytes_per_voxel)
    {
        case 1:
            if (is_signed)
                return CV_8SC(n_channels);
            return CV_8UC(n_channels);
        case 2:
            if (is_signed)
                return CV_16SC(n_channels);
            return CV_16UC(n_channels);
        case 4:
            if (is_floating_point)
                return CV_32FC(n_channels);
            if (!is_signed)
                MCP3D_MESSAGE("warning: opencv does not have CV_32U, returning CV_32S instead")
            return CV_32SC(n_channels);
        case 8:
            if (!is_floating_point)
                MCP3D_DOMAIN_ERROR("opencv does not have 64bit integer CvType")
            return CV_64FC(n_channels);
        default:
            MCP3D_INVALID_ARGUMENT("invalid bytes per voxel number")
    }
}

hid_t mcp3d::VoxelTypeToH5Type(VoxelType voxel_type)
{
    if (voxel_type == mcp3d::VoxelType::M8U)
        return H5T_NATIVE_UCHAR;
    else if (voxel_type == mcp3d::VoxelType::M8S)
        return H5T_NATIVE_CHAR;
    else if (voxel_type == mcp3d::VoxelType::M16U)
        return H5T_NATIVE_USHORT;
    else if (voxel_type == mcp3d::VoxelType::M16S)
        return H5T_NATIVE_SHORT;
    else if (voxel_type == mcp3d::VoxelType::M32U)
        return H5T_NATIVE_UINT32;
    else if (voxel_type == mcp3d::VoxelType::M32S)
        return H5T_NATIVE_INT32;
    else if (voxel_type == mcp3d::VoxelType::M64U)
        return H5T_NATIVE_UINT64;
    else if (voxel_type == mcp3d::VoxelType::M64S)
        return H5T_NATIVE_INT64;
    else if (voxel_type == mcp3d::VoxelType::M32F)
        return H5T_NATIVE_FLOAT;
    else if (voxel_type == mcp3d::VoxelType::M64F)
        return H5T_NATIVE_DOUBLE;
    else
        MCP3D_INVALID_ARGUMENT("voxel type can not be unknown")
}

