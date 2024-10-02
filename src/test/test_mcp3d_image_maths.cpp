//
// Created by muyezhu on 6/11/19.
//
#include <cstring>
#include <Eigen/Core>
#include <unsupported/Eigen/CXX11/Tensor>
#include <gtest/gtest.h>
#include "common/mcp3d_types.hpp"
#include "common/mcp3d_utility.hpp"
#include "image/mcp3d_image_maths.hpp"

using namespace std;

TEST(ImageMaths, CastVolume)
{
    vector<int> dims {100, 100, 3, 2};
    size_t n = mcp3d::ReduceProdSeq<size_t>(dims);
    // src_data holds uint16_t data
    unique_ptr<uint8_t []> src_data(new uint8_t[n * sizeof(uint16_t)]),
                           target_data;
    // random initialization of src_data
    mcp3d::MCPTensor1DMap<uint16_t> src_map((uint16_t*)(src_data.get()), n);
    src_map.setRandom();
    // make a copy of the src_data
    unique_ptr<uint8_t []> src_copy(new uint8_t[n * sizeof(uint16_t)]);
    memcpy(src_copy.get(), src_data.get(), n * sizeof(uint16_t));
    // create correct outputs into tensors
    mcp3d::MCPTensor1D<uint8_t> result_cast8 = src_map.cast<uint8_t>();
    mcp3d::MCPTensor1D<uint16_t> result_cast16 = src_map.cast<uint16_t>();
    mcp3d::MCPTensor1D<float> result_castf = src_map.cast<float>();
    // cast into target
    mcp3d::CastVolume<uint16_t, uint8_t>(src_data, dims, target_data);
    EXPECT_TRUE(memcmp(target_data.get(), result_cast8.data(), n * sizeof(uint8_t)) == 0);
    EXPECT_TRUE(memcmp(src_data.get(), src_copy.get(), n * sizeof(uint16_t)) == 0);
    mcp3d::CastVolume<uint16_t, uint16_t>(src_data, dims, target_data);
    EXPECT_TRUE(memcmp(target_data.get(), result_cast16.data(), n * sizeof(uint16_t)) == 0);
    EXPECT_TRUE(memcmp(src_data.get(), src_copy.get(), n * sizeof(uint16_t)) == 0);
    // cast in place
    mcp3d::CastVolume<uint16_t, uint16_t>(src_data, dims, src_data);
    EXPECT_TRUE(memcmp(src_data.get(), result_cast16.data(), n * sizeof(uint16_t)) == 0);
    EXPECT_TRUE(memcmp(src_data.get(), src_copy.get(), n * sizeof(uint16_t)) == 0);
    mcp3d::CastVolume<uint16_t, float>(src_data, dims, src_data);
    EXPECT_TRUE(memcmp(src_data.get(), result_castf.data(), n * sizeof(float)) == 0);
}

TEST(ImageMaths, NormalizeVolume)
{
    vector<int> dims {10, 1024, 1024};
    size_t n = mcp3d::ReduceProdSeq<size_t>(dims);
    unique_ptr<uint8_t []> src_data(new uint8_t[n * sizeof(uint16_t)]), target_data;
    // random initialization of uint16 src_data
    mcp3d::SetRandom<uint16_t >(src_data.get(), dims);
    // make a copy of the src_data
    unique_ptr<uint8_t []> src_copy(new uint8_t[n * sizeof(uint16_t)]);
    memcpy(src_copy.get(), src_data.get(), n * sizeof(uint16_t));
    // normalize to target
    mcp3d::NormalizeVolume<uint16_t, uint8_t>(src_data, dims, target_data);
    mcp3d::MCPTensor1DMap<uint8_t> target_map8(target_data.get(), n);
    mcp3d::MCPTensor0D<uint8_t> target_map8_min = target_map8.minimum();
    EXPECT_EQ(0, target_map8_min(0));
    mcp3d::MCPTensor0D<uint8_t> target_map8_max = target_map8.maximum();
    EXPECT_EQ(255, target_map8_max(0));
    EXPECT_TRUE(memcmp(src_data.get(), src_copy.get(), n * sizeof(uint16_t)) == 0);
    // normalize to src
    mcp3d::NormalizeVolume<uint16_t, float>(src_data, dims, src_data);
    mcp3d::MCPTensor1DMap<float> src_map_float((float*)(src_data.get()), n);
    mcp3d::MCPTensor0D<bool> in_bounds = (src_map_float >= 0.0f).all();
    EXPECT_TRUE(in_bounds(0));
    in_bounds = (src_map_float <= 1.0f).all();
    EXPECT_TRUE(in_bounds(0));
    // random initialization of double src_data
    src_data = make_unique<uint8_t []>(n * sizeof(double));
    mcp3d::MCPTensor1DMap<double> src_map_double((double*)(src_data.get()), n);
    src_map_double.setRandom();
    mcp3d::NormalizeVolume<double, double>(src_data, dims, src_data);
    in_bounds = (src_map_double >= 0.0).all();
    EXPECT_TRUE(in_bounds(0));
    in_bounds = (src_map_double <= 1.0).all();
    EXPECT_TRUE(in_bounds(0));
}

TEST(ImageMaths, TopPercentile)
{

}
