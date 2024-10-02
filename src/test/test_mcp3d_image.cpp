//
// Created by muyezhu on 2/28/18.
//
#include <gtest/gtest.h>
#include "image/mcp3d_image.hpp"
#include "image/mcp3d_image_in_memory.hpp"
#include "image_layout/mcp3d_image_layout_constants.hpp"

using namespace std;

string test_mimage_tiff_format_dir(const string& subdir = string())
{
    string d = mcp3d::JoinPath({mcp3d::test_data_dir(),
                                "image_formats", "tiff_format", subdir});
    MCP3D_ASSERT(mcp3d::IsDir(d))
    return d;
}

TEST(MImage, DefaultConstructor)
{
    mcp3d::MImage img;
    EXPECT_TRUE(img.can_read_storage());
    EXPECT_TRUE(img.image_info().empty());
    EXPECT_TRUE(img.channel_infos().empty());
    EXPECT_TRUE(img.selected_view().empty());
    EXPECT_TRUE(img.loaded_view().empty());
    EXPECT_EQ(vector<int>({0, 0, 0, 0, 0}), img.dims());
    EXPECT_EQ(vector<int>({0, 0, 0}), img.xyz_dims(0, 0));
    EXPECT_EQ(0, img.xdim());
    EXPECT_EQ(0, img.ydim());
    EXPECT_EQ(0, img.zdim());
    EXPECT_EQ(0, img.n_channels());
    EXPECT_EQ(0, img.n_times());
    EXPECT_EQ(0, img.n_volumes());
    EXPECT_EQ(mcp3d::VoxelType::UNKNOWN, img.voxel_type());
    EXPECT_TRUE(img.volume_root_dir().empty());
    EXPECT_EQ(0, img.n_pyr_levels());
    EXPECT_EQ(0, img.SelectedBytesPerVoxel());
    EXPECT_EQ(0, img.LoadedBytesPerVoxel());
    EXPECT_EQ(0, img.StorageBytesPerVoxel());
    EXPECT_FALSE(img.is_wrapper());
    EXPECT_THROW(img.channel_pyr_info(0, 0), mcp3d::MCPAssertionError);
}

TEST(MImageInMemory, DefaultConstructor)
{
    mcp3d::MImageInMemory img;
    EXPECT_FALSE(img.can_read_storage());
    EXPECT_TRUE(img.image_info().VolumeEmpty());
    EXPECT_TRUE(img.channel_infos().empty());
    EXPECT_TRUE(img.selected_view().empty());
    EXPECT_TRUE(img.loaded_view().empty());
    EXPECT_EQ(vector<int>({0, 0, 0, 0, 0}), img.dims());
    EXPECT_EQ(vector<int>({0, 0, 0}), img.xyz_dims(0, 0));
    EXPECT_EQ(0, img.xdim());
    EXPECT_EQ(0, img.ydim());
    EXPECT_EQ(0, img.zdim());
    EXPECT_EQ(0, img.n_channels());
    EXPECT_EQ(0, img.n_times());
    EXPECT_EQ(0, img.n_volumes());
    EXPECT_EQ(mcp3d::VoxelType::UNKNOWN, img.voxel_type());
    EXPECT_TRUE(img.image_root_dir().empty());
    EXPECT_EQ(0, img.n_pyr_levels());
    EXPECT_EQ(0, img.SelectedBytesPerVoxel());
    EXPECT_EQ(0, img.LoadedBytesPerVoxel());
    EXPECT_EQ(0, img.StorageBytesPerVoxel());
    EXPECT_FALSE(img.is_wrapper());
    EXPECT_THROW(img.channel_pyr_info(0, 0), mcp3d::MCPAssertionError);
}

TEST(MImageInMemory, ConstructFromDimensionType)
{
    vector<int> dims({789, 456, 123}), dims_5d = mcp3d::To5D(dims);
    mcp3d::VoxelType voxel_8u = mcp3d::VoxelType::M8U;
    uint8_t voxel_val0 = 127, voxel_val1 = 128;

    mcp3d::MImageInMemory img(dims, voxel_8u);
    EXPECT_EQ(dims_5d, img.image_info().dims());
    EXPECT_FALSE(img.selected_view().empty());
    EXPECT_FALSE(img.loaded_view().empty());
    EXPECT_EQ(1, img.loaded_view().n_volumes());
    EXPECT_EQ(img.selected_view(), img.loaded_view());
    EXPECT_EQ(dims, img.selected_view().xyz_dims(0, 0));
    EXPECT_EQ(voxel_8u, img.selected_view().voxel_type());
    EXPECT_FALSE(img.is_wrapper());
    EXPECT_EQ(mcp3d::VoxelSize(voxel_8u), img.LoadedBytesPerVoxel());
    EXPECT_EQ(mcp3d::ReduceProdSeq<int>(dims), img.loaded_view().VolumeBytes(0));

    uint8_t* ptr = img.Volume(0);
    ptr[0] = voxel_val0;
    EXPECT_EQ(voxel_val0, img.At(0, 0, 0, 0, 0));
    img.SetVoxel(0, 0, 0, 0, voxel_val1, 0);
    EXPECT_EQ(voxel_val1, img.At(0, 0, 0, 0, 0));
}

TEST(MImageInMemory, WrapData)
{
    vector<int> dims({3, 1024, 1024});
    long n_bytes_volume = mcp3d::ReduceProdSeq<size_t>(dims);
    unique_ptr<uint16_t []> data = make_unique<uint16_t []>((size_t)n_bytes_volume);
    mcp3d::SetRandom<uint16_t>(data.get(), dims);
    mcp3d::MImageInMemory wrapper;
    wrapper.WrapData(vector<uint16_t*>({data.get()}), dims);
    EXPECT_EQ(mcp3d::IN_MEMORY_STRING, wrapper.image_info().channel_root_dir());
    EXPECT_EQ(1, wrapper.n_pyr_levels());
    EXPECT_EQ(dims, wrapper.xyz_dims(0, 0));
    EXPECT_EQ(vector<int>({1, 1, 3, 1024, 1024}), wrapper.dims());
    EXPECT_THROW(wrapper.xyz_dims(1, 0), mcp3d::MCPAssertionError);
    EXPECT_THROW(wrapper.xyz_dims(0, 1), mcp3d::MCPAssertionError);
    EXPECT_FALSE(wrapper.can_read_storage());
    EXPECT_TRUE(wrapper.is_wrapper());
    EXPECT_FALSE(wrapper.selected_view().empty());
    EXPECT_EQ(wrapper.selected_view(), wrapper.loaded_view());
    EXPECT_EQ(mcp3d::VoxelType::M16U, wrapper.selected_view().voxel_type());
    EXPECT_EQ(dims, wrapper.selected_view().xyz_dims(0));
    EXPECT_EQ(2, wrapper.selected_view().VoxelBytes());
    EXPECT_EQ(n_bytes_volume, wrapper.selected_view().VolumeVoxels(0));
    for (long i = 0; i < n_bytes_volume; ++i)
        EXPECT_EQ(data.get()[i], wrapper.Volume<uint16_t>()[i]);
}

TEST(MImage, MaxPyramidLevel)
{
    mcp3d::MImage img1(test_mimage_tiff_format_dir("large_files_read"));
    img1.ReadImageInfo(0);
    EXPECT_EQ(7, img1.MaxPyramidLevel(0));
    EXPECT_EQ(8, img1.MaxNumPyramidLevels(0));

    vector<int> dims({10, 10, 10});
    mcp3d::VoxelType voxel_8u = mcp3d::VoxelType::M8U;
    mcp3d::MImageInMemory img2(dims, voxel_8u);
    EXPECT_EQ(0, img2.MaxPyramidLevel(0));
    EXPECT_EQ(1, img2.MaxNumPyramidLevels(0));
}

TEST(MImage, PyramidBoundary)
{
    mcp3d::MImage img(test_mimage_tiff_format_dir("large_files_read"));
    // x dimension 22343, y dimension 26604
    img.ReadImageInfo(0);
    vector<int> xyz_dims0 = img.xyz_dims(0, 0);
    vector<int> xyz_dims1 = img.xyz_dims(0, 1);
    mcp3d::MImageBlock block_level0{vector<int>({0, xyz_dims0[1] - 1, xyz_dims0[2] - 1}),
                                    vector<int>({xyz_dims0[0] + 1, 100, 10})};
    mcp3d::MImageBlock block_level1{vector<int>({0, xyz_dims1[1] - 1, xyz_dims1[2]}),
                                    vector<int>({xyz_dims1[0] + 1, 50, 50})};
    EXPECT_NO_THROW(img.SelectView(block_level0, 0, 0, false));
    // global block offsets are within boundary
    EXPECT_NO_THROW(img.SelectView(block_level0, 0, 1, false));
    EXPECT_NO_THROW(img.SelectView(block_level1, 0, 1, false));
    EXPECT_NO_THROW(img.SelectView(block_level1, 0, 1, true));

    img.SelectView(block_level1, 0, 1, true);
    EXPECT_TRUE(img.selected_view().OutOfPyrImageBoundary(img.image_info()));
    // TODO: construct cases thats only partially out of bounds
    // change expect_false below to expect_true: fully out of bounds no longer returns
    // false in PartiallyOutOfPyrImageBoundary
    EXPECT_FALSE(img.selected_view().PartiallyOutOfPyrImageBoundary(img.image_info()));
    img.ReadData();
    for (long i = 0; i < img.loaded_view().VolumeVoxels(); ++i)
        EXPECT_EQ(0, img.Volume()[i]);
}

TEST(MImage, ReadImageInfoFromEmptyDir)
{
    string empty_dir(mcp3d::JoinPath({mcp3d::test_data_dir(), "MImage_ReadImageInfoFromEmptyDir"}));
    EXPECT_TRUE(mcp3d::RemovePath(empty_dir));
    mcp3d::MakeDirectories(empty_dir);
    mcp3d::MImage img(empty_dir);
    EXPECT_NO_THROW(img.ReadImageInfo(0));
    EXPECT_FALSE(img.image_info().empty());
    EXPECT_EQ(empty_dir, img.image_info().volume_root_dir());
    EXPECT_FALSE(img.channel_info(0).empty());
    EXPECT_EQ(empty_dir, img.channel_info(0).channel_root_dir());
    EXPECT_TRUE(img.image_info().channel_info_read(0));
}

