//
// Created by muyezhu on 2/12/18.
//
#include <unordered_map>
#include <boost/filesystem.hpp>
#include <gtest/gtest.h>
#include "common/mcp3d_common.hpp"
#include "project_structure/mcp3d_project_structure_constants.hpp"
#include "project_structure/mcp3d_channel_layout.hpp"
#include "image/mcp3d_channel_info.hpp"
#include "image/mcp3d_image.hpp"

using namespace std;

string test_img_info_dir()
{
    string result(mcp3d::JoinPath({mcp3d::test_data_dir(), "image_info"}));
    MCP3D_ASSERT(mcp3d::IsDir(result))
    return result;
}

string test_channel_pyr_info_json_path()
{
    string result(mcp3d::JoinPath({mcp3d::test_data_dir(), "image_info", "test_channel_pyr_info.json"}));
    MCP3D_ASSERT(mcp3d::IsFile(result))
    return result;
}

string test_channel_info_json_path()
{
    string result(mcp3d::JoinPath({mcp3d::test_data_dir(), "image_info", "test_channel_info.json"}));
    MCP3D_ASSERT(mcp3d::IsFile(result))
    return result;
}

string img_pyramids_dir()
{
    string result = mcp3d::JoinPath({mcp3d::test_data_dir(), "image_formats", "tiff_format", "large_files_read"});
    MCP3D_ASSERT(mcp3d::IsDir(result))
    return result;
}

TEST(MChannelPyrInfo, DefaultConstructor)
{
    mcp3d::MChannelPyrInfo channel_pyr_info {};
    EXPECT_TRUE(channel_pyr_info.channel_pyr_slices().IsFlat());
    EXPECT_EQ(mcp3d::FileFormat::UNKNOWN, channel_pyr_info.file_format());
    EXPECT_EQ(mcp3d::VoxelType::UNKNOWN, channel_pyr_info.voxel_type());
    EXPECT_FALSE(channel_pyr_info.is_level0());
    EXPECT_FALSE(channel_pyr_info.zdim_scaled());
    EXPECT_TRUE(channel_pyr_info.channel_pyr_dir().empty());
    EXPECT_TRUE(channel_pyr_info.slice_image_names().empty());
    EXPECT_TRUE(channel_pyr_info.dims_order().empty());
    EXPECT_EQ(0, channel_pyr_info.xdim());
    EXPECT_EQ(0, channel_pyr_info.ydim());
    EXPECT_EQ(0, channel_pyr_info.zdim());
    EXPECT_EQ(0, channel_pyr_info.chunk_xdim());
    EXPECT_EQ(0, channel_pyr_info.chunk_ydim());
    EXPECT_EQ(0, channel_pyr_info.chunk_zdim());
    EXPECT_EQ(0, channel_pyr_info.n_xchunks());
    EXPECT_EQ(0, channel_pyr_info.n_ychunks());
    EXPECT_EQ(0, channel_pyr_info.n_zchunks());
    EXPECT_EQ(0, channel_pyr_info.n_total_chunks());
}

TEST(MChannelPyrInfo, CopyConstruct)
{
    mcp3d::MChannelPyrInfo channel_pyr_info0(test_channel_pyr_info_json_path());
    mcp3d::MChannelPyrInfo copy0(channel_pyr_info0);
    EXPECT_TRUE(channel_pyr_info0 == copy0);
    mcp3d::MChannelInfo channel_info(channel_pyr_info0);
    const mcp3d::MChannelPyrInfo& channel_pyr_info1 = channel_info.channel_pyr_infos()[0];
    mcp3d::MChannelPyrInfo copy1(channel_pyr_info1);
    EXPECT_TRUE(channel_pyr_info1 == copy1);
}

TEST(MChannelPyrInfo, ConstructFromFile)
{
    mcp3d::MChannelPyrInfo channel_pyr_info(test_channel_pyr_info_json_path());
    EXPECT_EQ(mcp3d::FileFormat::TIFF, channel_pyr_info.file_format());
    EXPECT_EQ(mcp3d::VoxelType::M8U, channel_pyr_info.voxel_type());
    EXPECT_TRUE(size_t(channel_pyr_info.n_total_chunks()) ==
                        channel_pyr_info.slice_image_names().at("").size());
    EXPECT_EQ(mcp3d::ParentDir(test_channel_pyr_info_json_path()),
              channel_pyr_info.channel_pyr_dir());
    vector<string> tiff_names({"a_x0_y0_z0.tif", "a_x1024_y0_z0.tif",
                               "a_x0_y1024_z0.tif", "a_x1024_y1024_z0.tif"});
    unordered_map<string, vector<string>> slice_tiff_names({{"", tiff_names}});
    EXPECT_EQ(slice_tiff_names, channel_pyr_info.slice_image_names());
    EXPECT_EQ(vector<int>({1024, 2048, 2048}), channel_pyr_info.xyz_dims());
    EXPECT_EQ(2048, channel_pyr_info.xdim());
    EXPECT_EQ(2048, channel_pyr_info.ydim());
    EXPECT_EQ(1024, channel_pyr_info.zdim());
    EXPECT_EQ(2048, channel_pyr_info.slice_dim(mcp3d::ChannelAxis::X));
    EXPECT_EQ(2048, channel_pyr_info.slice_dim(mcp3d::ChannelAxis::Y));
    EXPECT_EQ(1024, channel_pyr_info.slice_dim(mcp3d::ChannelAxis::Z));
    EXPECT_EQ(1, channel_pyr_info.slice_number(mcp3d::ChannelAxis::X));
    EXPECT_EQ(1, channel_pyr_info.slice_number(mcp3d::ChannelAxis::Y));
    EXPECT_EQ(1, channel_pyr_info.slice_number(mcp3d::ChannelAxis::Z));
    EXPECT_EQ(1024, channel_pyr_info.chunk_xdim());
    EXPECT_EQ(1024, channel_pyr_info.chunk_ydim());
    EXPECT_EQ(1024, channel_pyr_info.chunk_zdim());
    EXPECT_EQ(2, channel_pyr_info.n_xchunks());
    EXPECT_EQ(2, channel_pyr_info.n_ychunks());
    EXPECT_EQ(1, channel_pyr_info.n_zchunks());
    EXPECT_EQ(4, channel_pyr_info.n_total_chunks());
    EXPECT_EQ("zyx", channel_pyr_info.dims_order());
    EXPECT_TRUE(channel_pyr_info.is_level0());
    EXPECT_FALSE(channel_pyr_info.zdim_scaled());
}

TEST(MChannelInfo, DefaultConstructor)
{
    mcp3d::MChannelInfo info;
    EXPECT_TRUE(info.channel_root_dir().empty());
    EXPECT_EQ(mcp3d::ZSCALE_NONE, info.z_scale_start_level());
    EXPECT_TRUE(info.channel_pyr_infos().empty());
    EXPECT_TRUE(info.pyr_xy_ratios().empty());
    EXPECT_TRUE(info.pyr_z_ratios().empty());
}

TEST(MChannelInfo, ConstructFromDimensionAndVoxelType)
{
    vector<int> dimensions({10, 1024, 1024});
    mcp3d::MChannelInfo info_from_memory(dimensions, mcp3d::VoxelType::M32F);
    EXPECT_EQ(mcp3d::IN_MEMORY_STRING, info_from_memory.channel_root_dir());
    EXPECT_EQ((int)sizeof(float), info_from_memory.VoxelBytes(0));
    EXPECT_EQ(mcp3d::VoxelType::M32F, info_from_memory.voxel_type(0));
    EXPECT_FALSE(info_from_memory.empty());
    EXPECT_EQ(1, info_from_memory.n_pyr_infos());
    EXPECT_THROW(info_from_memory.ImagePath(1, 0), mcp3d::MCPAssertionError);
    EXPECT_EQ(mcp3d::IN_MEMORY_STRING, info_from_memory.ImagePath(0, 1));
    EXPECT_EQ(mcp3d::IN_MEMORY_STRING, info_from_memory.ImagePath(0, 0));
    EXPECT_THROW(info_from_memory.xyz_dims(1), mcp3d::MCPAssertionError);
    EXPECT_THROW(info_from_memory.xyz_chunk_dims(1), mcp3d::MCPAssertionError);
    EXPECT_EQ(dimensions, info_from_memory.xyz_dims(0));
    EXPECT_EQ(dimensions, info_from_memory.xyz_chunk_dims(0));
    EXPECT_EQ(vector<int>({1}), info_from_memory.pyr_xy_ratios());
    EXPECT_EQ((size_t)1, info_from_memory.channel_pyr_infos().size());
    EXPECT_EQ(mcp3d::IN_MEMORY_STRING, info_from_memory.channel_pyr_info(0).channel_pyr_dir());
    EXPECT_EQ(1, info_from_memory.channel_pyr_info(0).n_total_chunks());
    EXPECT_EQ("zyx", info_from_memory.channel_pyr_info(0).dims_order());
    EXPECT_TRUE(info_from_memory.channel_pyr_info(0).is_level0());
    EXPECT_EQ(mcp3d::ZSCALE_NONE, info_from_memory.z_scale_start_level());
    EXPECT_EQ(mcp3d::FileFormat::UNKNOWN, info_from_memory.file_format(0));
}

TEST(MChannelInfo, Clear)
{
    vector<int> dimensions({10, 1024, 1024});
    mcp3d::MChannelInfo info_from_memory(dimensions, mcp3d::VoxelType::M32F);
    info_from_memory.Clear();
    EXPECT_EQ(mcp3d::MChannelInfo{}, info_from_memory);
}

TEST(MChannelInfo, ConstructFromPyrInfos)
{
    vector<int> dimensions({10, 1024, 1024});
    mcp3d::MChannelInfo from_memory_info(dimensions, mcp3d::VoxelType::M32F);
    vector<mcp3d::MChannelPyrInfo> pyr_infos;
    pyr_infos.push_back(from_memory_info.channel_pyr_infos()[0]);
    EXPECT_THROW(mcp3d::MChannelInfo{pyr_infos}, mcp3d::MCPAssertionError);

    mcp3d::MChannelPyrInfo pyr_info0{mcp3d::JoinPath({test_img_info_dir(), "same_root", "__pyr_level_00_info__.json"})},
                           pyr_info1{mcp3d::JoinPath({test_img_info_dir(), "same_root", "pyr_level_01", "__pyr_level_01_info__.json"})},
                           pyr_info2{mcp3d::JoinPath({test_img_info_dir(), "same_root", "__pyr_level_02_info__.json"})};

    pyr_infos.clear();
    pyr_infos.emplace_back(mcp3d::JoinPath({test_img_info_dir(), "same_root", "__pyr_level_00_info__.json"}));
    pyr_infos.emplace_back(mcp3d::JoinPath({test_img_info_dir(), "same_root", "__pyr_level_02_info__.json"}));
    pyr_infos.emplace_back(mcp3d::JoinPath({test_img_info_dir(), "same_root", "pyr_level_01", "__pyr_level_01_info__.json"}));
    mcp3d::MChannelInfo same_root_info(pyr_infos);
    EXPECT_EQ(mcp3d::JoinPath({test_img_info_dir(), "same_root"}),
              same_root_info.channel_root_dir());
    EXPECT_EQ(mcp3d::JoinPath({test_img_info_dir(), "same_root"}),
              same_root_info.channel_pyr_info(0).channel_pyr_dir());
    EXPECT_EQ(mcp3d::JoinPath({test_img_info_dir(), "same_root", "pyr_level_01"}),
              same_root_info.channel_pyr_info(1).channel_pyr_dir());
    EXPECT_EQ(mcp3d::JoinPath({test_img_info_dir(), "same_root"}),
              same_root_info.channel_pyr_info(2).channel_pyr_dir());
    EXPECT_EQ(vector<int>({1, 2, 4}), same_root_info.pyr_xy_ratios());
    EXPECT_EQ(vector<int>({1, 1, 2}), same_root_info.pyr_z_ratios());
    EXPECT_EQ(2, same_root_info.z_scale_start_level());
    EXPECT_EQ(mcp3d::VoxelType::M16U, same_root_info.voxel_type(1));
    EXPECT_EQ(vector<int>({3, 6651, 5585}), same_root_info.xyz_dims(2));
    EXPECT_EQ(vector<int>({1, 6651, 5585}), same_root_info.xyz_chunk_dims(2));

    pyr_infos.clear();
    pyr_infos.emplace_back(mcp3d::JoinPath({test_img_info_dir(), "not_same_root0", "level0", "__pyr_level_00_info__.json"}));
    pyr_infos.emplace_back(mcp3d::JoinPath({test_img_info_dir(), "not_same_root1", "1", "__pyr_level_01_info__.json"}));
    pyr_infos.emplace_back(mcp3d::JoinPath({test_img_info_dir(), "not_same_root1", "2", "__pyr_level_02_info__.json"}));
    mcp3d::MChannelInfo not_same_root_info(pyr_infos);
    EXPECT_EQ(mcp3d::JoinPath({test_img_info_dir(), "not_same_root0", "level0"}),
              not_same_root_info.channel_root_dir());
    EXPECT_EQ(mcp3d::JoinPath({test_img_info_dir(), "not_same_root0", "level0"}),
              not_same_root_info.channel_pyr_infos()[0].channel_pyr_dir());
    EXPECT_EQ(mcp3d::JoinPath({test_img_info_dir(), "not_same_root1", "1"}),
              not_same_root_info.channel_pyr_infos()[1].channel_pyr_dir());
    EXPECT_EQ(mcp3d::JoinPath({test_img_info_dir(), "not_same_root1", "2"}),
              not_same_root_info.channel_pyr_infos()[2].channel_pyr_dir());
    EXPECT_EQ(pyr_info1.xyz_dims(), not_same_root_info.xyz_dims(1));
    EXPECT_EQ(vector<int>({3, 6651, 5585}), not_same_root_info.xyz_dims(2));
    EXPECT_EQ(vector<int>({1, 6651, 5585}), not_same_root_info.xyz_chunk_dims(2));
    EXPECT_EQ(vector<int>({1, 2, 4}), not_same_root_info.pyr_xy_ratios());
    EXPECT_EQ(vector<int>({1, 1, 2}), not_same_root_info.pyr_z_ratios());
    EXPECT_EQ(2, not_same_root_info.z_scale_start_level());
    EXPECT_EQ(mcp3d::VoxelType::M16U, not_same_root_info.voxel_type(1));
}

TEST(MChannelInfo, SerializeAndDeserialize)
{
    string channel_info_path(test_channel_info_json_path());
    mcp3d::MChannelInfo channel_info(channel_info_path);
    channel_info.Save();
    // read -> save -> read should produce equal MChannelInfo objects
    mcp3d::MChannelInfo channel_info_save(channel_info_path);
    EXPECT_TRUE(channel_info == channel_info_save);
}

TEST(MChannelInfo, SaveToFile)
{
    mcp3d::MImage image(img_pyramids_dir());
    string channel_info_path (mcp3d::MChannelLayout::channel_info_path(img_pyramids_dir(), 0));
    mcp3d::RemovePath(channel_info_path);
    image.ReadImageInfo(0);
    image.SaveImageInfo();

    // written image info json file exists
    EXPECT_TRUE(mcp3d::IsFile(channel_info_path));
    // read -> save -> read should produce equal MChannelInfo objects
    mcp3d::MChannelInfo saved_channel_info(channel_info_path);
    EXPECT_TRUE(image.channel_info(0) == saved_channel_info);
}



