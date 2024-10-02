//
// Created by muyezhu on 1/16/19.
//
#include <gtest/gtest.h>
#include "common/mcp3d_exceptions.hpp"
#include "image/mcp3d_image_info.hpp"
#include "image/mcp3d_image_view.hpp"

using namespace std;

string test_image_view_dir()
{
    string d(mcp3d::JoinPath({mcp3d::test_data_dir(), "image_view"}));
    MCP3D_ASSERT(mcp3d::IsDir(d))
    return d;
}

TEST(MImageBlock, Constructor)
{
    // invalid constructor arguments
    EXPECT_THROW(mcp3d::MImageBlock(vector<int>({0, 0, 0, 0})),
                 mcp3d::MCPAssertionError);
    EXPECT_THROW(mcp3d::MImageBlock(vector<int>({0, 0, -1})),
                 mcp3d::MCPOutOfRangeError);
    EXPECT_THROW(mcp3d::MImageBlock(vector<int>({0, 0, 0}), vector<int>({0, 0, -1})),
                 mcp3d::MCPOutOfRangeError);
    EXPECT_THROW(mcp3d::MImageBlock(vector<int>({0, 0, 0}), vector<int>({0, 0, 0}), vector<int>({1, 1, 0})),
                 mcp3d::MCPOutOfRangeError);

    // default constructor
    mcp3d::MImageBlock block_default {};
    EXPECT_EQ(vector<int>({0, 0, 0}), block_default.offsets());
    EXPECT_EQ(vector<int>({0, 0, 0}), block_default.extents());
    EXPECT_EQ(vector<int>({1, 1, 1}), block_default.strides());

    // general constructors
    mcp3d::MImageBlock block0(vector<int>({100, 100}), vector<int>({0, 10}), vector<int>({2}));
    EXPECT_EQ(vector<int>({0, 100, 100}), block0.offsets());
    EXPECT_EQ(vector<int>({0, 0, 10}), block0.extents());
    EXPECT_EQ(vector<int>({1, 1, 2}), block0.strides());
}

TEST(MImageBlock, ClearBlock)
{
    mcp3d::MImageBlock block(vector<int>({100, 100}), vector<int>({0, 10}), vector<int>({2}));
    block.Clear();
    EXPECT_EQ(vector<int>({0, 0, 0}), block.offsets());
    EXPECT_EQ(vector<int>({0, 0, 0}), block.extents());
    EXPECT_EQ(vector<int>({1, 1, 1}), block.strides());
}

TEST(MImageBlock, Equality)
{
    mcp3d::MImageBlock block0(vector<int>({100, 100}), vector<int>({0, 10}), vector<int>({2}));
    mcp3d::MImageBlock block1(vector<int>({0, 100, 100}), vector<int>({0, 10}), vector<int>({1, 1, 2}));
    mcp3d::MImageBlock block2(vector<int>({0, 100, 100}), vector<int>({0, 10}), vector<int>({1, 2, 2}));
    EXPECT_TRUE(block0 == block1);
    EXPECT_FALSE(block1 == block2);
    block0.Clear();
    EXPECT_TRUE(block0 == mcp3d::MImageBlock{});
}

TEST(MImageBlock, Copy)
{
    mcp3d::MImageBlock block0(vector<int>({100, 100}), vector<int>({0, 10}), vector<int>({2}));
    mcp3d::MImageBlock block1 {};
    EXPECT_FALSE(block0 == block1);
    mcp3d::MImageBlock block2 = block0;
    EXPECT_TRUE(block0 == block2);
    block1 = block2;
    EXPECT_TRUE(block0 == block1);
}

TEST(MImageView, Constructor)
{
    mcp3d::MImageView view {};
    EXPECT_TRUE(view.global_image_block().empty());
    EXPECT_TRUE(view.empty());
    EXPECT_EQ(mcp3d::VoxelType::UNKNOWN, view.voxel_type());
    EXPECT_EQ(-1, view.pyr_level());
    EXPECT_FALSE(view.interpret_block_as_local());
    EXPECT_EQ(vector<int>(3, -1), view.pyr_level_offsets(0));
    EXPECT_EQ(vector<int>(3, 0), view.pyr_level_extents(0));
    EXPECT_EQ(vector<int>(3, 0), view.pyr_level_strides(0));
}

TEST(MImageView, SelectView)
{
    string channel_info_path(mcp3d::JoinPath({test_image_view_dir(),
                                             "imaris_image_info.json"}));
    mcp3d::MChannelInfo channel_info(channel_info_path);
    channel_info.set_channel_number(0);
    mcp3d::MImageInfo image_info(channel_info.channel_root_dir(), {""});
    image_info += channel_info;
    mcp3d::MImageView image_view;
    mcp3d::MImageBlock block_all;
    image_view.SelectView(image_info, block_all, channel_info.channel_number(), 0, false);
    EXPECT_EQ(vector<int>({292, 6197, 5183}), image_view.xyz_dims(0));
    image_view.SelectView(image_info, block_all, channel_info.channel_number(), 1, false);
    EXPECT_EQ(vector<int>({292, 3098, 2591}), image_view.xyz_dims(0));
    image_view.SelectView(image_info, block_all, channel_info.channel_number(), 2, false);
    EXPECT_EQ(vector<int>({146, 1549, 1295}), image_view.xyz_dims(0));

    mcp3d::MImageBlock block0(vector<int>({100, 2000, 2000}), vector<int>({64, 2048, 2048}));
    image_view.SelectView(image_info, block0, channel_info.channel_number(), 0, false);
    EXPECT_EQ(vector<int>({100, 2000, 2000}), image_view.pyr_level_offsets(channel_info.channel_number()));
    image_view.SelectView(image_info, block0, channel_info.channel_number(), 1, false);
    EXPECT_EQ(vector<int>({100, 1000, 1000}), image_view.pyr_level_offsets(channel_info.channel_number()));
    image_view.SelectView(image_info, block0, channel_info.channel_number(), 2, false);
    EXPECT_EQ(vector<int>({50, 500, 500}), image_view.pyr_level_offsets(channel_info.channel_number()));
    image_view.SelectView(image_info, block0, channel_info.channel_number(), 3, false);
    EXPECT_EQ(vector<int>({25, 250, 250}), image_view.pyr_level_offsets(channel_info.channel_number()));

    mcp3d::MImageBlock block1(vector<int>({100, 2000, 2000}),
                              vector<int>({64, 2048, 2048}),
                              vector<int>({8, 8, 8}));
    image_view.SelectView(image_info, block1, channel_info.channel_number(), 0, false);
    EXPECT_EQ(vector<int>({100, 2000, 2000}), image_view.pyr_level_offsets(
            channel_info.channel_number()));
    EXPECT_EQ(vector<int>({64, 2048, 2048}), image_view.pyr_level_extents(
            channel_info.channel_number()));
    EXPECT_EQ(vector<int>({8, 8, 8}), image_view.pyr_level_strides(
            channel_info.channel_number()));
    image_view.SelectView(image_info, block1, channel_info.channel_number(), 1, false);
    EXPECT_EQ(vector<int>({100, 1000, 1000}), image_view.pyr_level_offsets(
            channel_info.channel_number()));
    EXPECT_EQ(vector<int>({64, 1024, 1024}), image_view.pyr_level_extents(
            channel_info.channel_number()));
    EXPECT_EQ(vector<int>({8, 4, 4}), image_view.pyr_level_strides(
            channel_info.channel_number()));
    image_view.SelectView(image_info, block1, channel_info.channel_number(), 2, false);
    EXPECT_EQ(vector<int>({50, 500, 500}), image_view.pyr_level_offsets(
            channel_info.channel_number()));
    EXPECT_EQ(vector<int>({32, 512, 512}), image_view.pyr_level_extents(
            channel_info.channel_number()));
    EXPECT_EQ(vector<int>({4, 2, 2}), image_view.pyr_level_strides(
            channel_info.channel_number()));
    image_view.SelectView(image_info, block1, channel_info.channel_number(), 3, false);
    EXPECT_EQ(vector<int>({25, 250, 250}), image_view.pyr_level_offsets(
            channel_info.channel_number()));
    EXPECT_EQ(vector<int>({16, 256, 256}), image_view.pyr_level_extents(
            channel_info.channel_number()));
    EXPECT_EQ(vector<int>({2, 1, 1}), image_view.pyr_level_strides(
            channel_info.channel_number()));
    image_view.SelectView(image_info, block1, channel_info.channel_number(), 4, false);
    EXPECT_EQ(vector<int>({12, 125, 125}), image_view.pyr_level_offsets(
            channel_info.channel_number()));
    EXPECT_EQ(vector<int>({8, 128, 128}), image_view.pyr_level_extents(
            channel_info.channel_number()));
    EXPECT_EQ(vector<int>({1, 1, 1}), image_view.pyr_level_strides(
            channel_info.channel_number()));

    mcp3d::MImageBlock block2(vector<int>({10, 200, 200}),
                              vector<int>({20, 200, 200}),
                              vector<int>({2, 2, 2}));
    image_view.SelectView(image_info, block2, channel_info.channel_number(), 1, true);
    EXPECT_EQ(vector<int>({10, 200, 200}), image_view.pyr_level_offsets(
            channel_info.channel_number()));
    EXPECT_EQ(vector<int>({20, 200, 200}), image_view.pyr_level_extents(
            channel_info.channel_number()));
    EXPECT_EQ(vector<int>({2, 2, 2}), image_view.pyr_level_strides(
            channel_info.channel_number()));
    EXPECT_EQ(vector<int>({10, 400, 400}), image_view.global_image_block().offsets());
    EXPECT_EQ(vector<int>({20, 400, 400}), image_view.global_image_block().extents());
    EXPECT_EQ(vector<int>({2, 4, 4}), image_view.global_image_block().strides());
    image_view.SelectView(image_info, block2, channel_info.channel_number(), 2, true);
    EXPECT_EQ(vector<int>({10, 200, 200}), image_view.pyr_level_offsets(
            channel_info.channel_number()));
    EXPECT_EQ(vector<int>({20, 200, 200}), image_view.pyr_level_extents(
            channel_info.channel_number()));
    EXPECT_EQ(vector<int>({2, 2, 2}), image_view.pyr_level_strides(
            channel_info.channel_number()));
    EXPECT_EQ(vector<int>({20, 800, 800}), image_view.global_image_block().offsets());
    EXPECT_EQ(vector<int>({40, 800, 800}), image_view.global_image_block().extents());
    EXPECT_EQ(vector<int>({4, 8, 8}), image_view.global_image_block().strides());

}

TEST(MImageView, OutOfBoundary)
{
    vector<mcp3d::MChannelPyrInfo> pyr_infos;
    // level 0 zyx dims 6, 26604, 22343
    // level 1 zyx dims 6, 13302, 11171
    pyr_infos.emplace_back(mcp3d::JoinPath({test_image_view_dir(),
                                            "__pyr_level_00_info__.json"}));
    pyr_infos.emplace_back(mcp3d::JoinPath({test_image_view_dir(),
                                            "__pyr_level_01_info__.json"}));

    mcp3d::MChannelInfo channel_info(pyr_infos);
    mcp3d::MImageInfo image_info(channel_info.channel_root_dir(), {""});
    image_info += channel_info;
    mcp3d::MImageView view {};

    // entirely out of bounds on level 1, legal on global level
    mcp3d::MImageBlock block0(vector<int>({0, 0, 11171}), vector<int>({100, 100, 100}));
    view.SelectView(image_info, block0, 0, 1, true);
    EXPECT_TRUE(view.OutOfPyrImageBoundary(image_info));
    EXPECT_TRUE(view.PartiallyOutOfPyrImageBoundary(image_info));

    // partially out of bounds
    mcp3d::MImageBlock block1(vector<int>({0, 20000, 20000}), vector<int>({10000, 10000, 10000}));
    view.SelectView(image_info, block1, 0, 0);
    EXPECT_FALSE(view.OutOfPyrImageBoundary(image_info));
    EXPECT_TRUE(view.PartiallyOutOfPyrImageBoundary(image_info));

    // local level block
    // offset along x allows a single voxel to be read from level 1.
    // the extent along x is 2. given a stride of 2, the view is in bound
    mcp3d::MImageBlock block2(vector<int>({0, 0, 11170}), vector<int>({1, 100, 2}), vector<int>({1, 1, 2}));
    view.SelectView(image_info, block2, 0, 1, true);
    EXPECT_FALSE(view.OutOfPyrImageBoundary(image_info));
    EXPECT_FALSE(view.PartiallyOutOfPyrImageBoundary(image_info));
}

