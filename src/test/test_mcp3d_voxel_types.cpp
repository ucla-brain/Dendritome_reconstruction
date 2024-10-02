//
// Created by muyezhu on 7/29/19.
//

#include <gtest/gtest.h>
#include "image_interface/mcp3d_voxel_types.hpp"

using namespace std;

TEST(VoxelTypes, VoxelTypeZeroValue)
{
    EXPECT_EQ(uint8_t(0) , mcp3d::VoxelTypeZeroValue<uint8_t>());
    EXPECT_EQ(int8_t(0), mcp3d::VoxelTypeZeroValue<int8_t>());
    EXPECT_EQ(uint16_t(0), mcp3d::VoxelTypeZeroValue<uint16_t>());
    EXPECT_EQ(int16_t(0), mcp3d::VoxelTypeZeroValue<int16_t>());
    EXPECT_EQ(uint32_t(0), mcp3d::VoxelTypeZeroValue<uint32_t>());
    EXPECT_EQ(int32_t(0), mcp3d::VoxelTypeZeroValue<int32_t>());
    EXPECT_EQ(uint64_t(0), mcp3d::VoxelTypeZeroValue<uint64_t>());
    EXPECT_EQ(int64_t(0), mcp3d::VoxelTypeZeroValue<int64_t>());
    EXPECT_EQ(float(0), mcp3d::VoxelTypeZeroValue<float_t>());
    EXPECT_EQ(double(0), mcp3d::VoxelTypeZeroValue<double_t>());
}

