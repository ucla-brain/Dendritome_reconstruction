//
// Created by muyezhu on 10/2/17.
//

#ifndef MCP3D_PIXEL_CONVERSION_HPP
#define MCP3D_PIXEL_CONVERSION_HPP

#include "common/mcp3d_common.hpp"
#include "mcp3d_image_legacy.hpp"

namespace mcp3d
{
class Image3D;

template <typename T>
void saturation_cast(const MCPTensor3D<T>& input,
                     Image3D& output);
}

template <typename T>
void mcp3d::saturation_cast(const mcp3d::MCPTensor3D<T>& input,
                            mcp3d::Image3D& output)
{
    if (! output.SameDimensions(input)) MCP3D_INVALID_ARGUMENT(
            "output must have same dimensions as input")

    T max_val;
    if (output.datatype() == VoxelType::M8U)
    {
        max_val = UCHAR_MAX;
        mcp3d::MCPTensor3D<T> saturate = input.cwiseMin(max_val);
        saturate = saturate.round();
        output.image<uint8_t>() = saturate.template cast<uint8_t>();
    }
    else if (output.datatype() == VoxelType::M16U)
    {
        max_val = USHRT_MAX;
        mcp3d::MCPTensor3D<T> saturate = input.cwiseMin(max_val);
        saturate = saturate.round();
        output.image<uint16_t>() = saturate.template cast<uint16_t>();
    }
    else if (output.datatype() == VoxelType::M32U)
    {
        max_val = ULONG_MAX;
        mcp3d::MCPTensor3D<T> saturate = input.cwiseMin(max_val);
        saturate = saturate.round();
        output.image<uint32_t>() = saturate.template cast<uint32_t>();
    }
    else if (output.datatype() == VoxelType::M32S)
    {
        max_val = LONG_MAX;
        mcp3d::MCPTensor3D<T> saturate = input.cwiseMin(max_val);
        saturate = saturate.round();
        output.image<int>() = saturate.template cast<int>();
    }
    else if (output.datatype() == VoxelType::M32F)
    {
        max_val = FLT_MAX;
        mcp3d::MCPTensor3D<T> saturate = input.cwiseMin(max_val);
        output.image<float>() = saturate.template cast<float>();
    }
    else if (output.datatype() == VoxelType::M64F)
    {
        max_val = DBL_MAX;
        mcp3d::MCPTensor3D<T> saturate = input.cwiseMin(max_val);
        output.image<double>() = saturate.template cast<double>();
    }
    else MCP3D_DOMAIN_ERROR("unknown data type")
}

#endif //MCP3D_PIXEL_CONVERSION_HPP
