//
// Created by muyezhu on 7/24/19.
//
#include <random>
#include <chrono>
#include <gtest/gtest.h>
#include "image/mcp3d_channel_info.hpp"
#include "image/mcp3d_image.hpp"
#include "image/mcp3d_ome_tiff_format.hpp"

using namespace std;

string test_hdf5_format_dir(const string& subdir = string());

string test_ome_tiff_format_dir(const string& subdir = string())
{
    string d = mcp3d::JoinPath({mcp3d::test_data_dir(),
                                "image_formats", "ome_tiff_format", subdir});
    MCP3D_ASSERT(mcp3d::IsDir(d))
    return d;
}

TEST(MOmeTiffFormat, ReadChannelInfoPyrDir)
{
    mcp3d::MOmeTiffFormat ome_format;
    // flat. data same as test_data/image_formats/tiff_format/Z0.tif, Z1.tif, Z2.tif
    string flat_dir(test_ome_tiff_format_dir("flat"));
    mcp3d::MChannelInfo channel_info = ome_format.ReadChannelInfo(flat_dir,
                                                                  0, true);
    EXPECT_EQ(flat_dir, channel_info.channel_root_dir());
    EXPECT_EQ(flat_dir, channel_info.channel_pyr_dir(0));
    EXPECT_EQ(1, channel_info.n_pyr_infos());
    EXPECT_EQ(2048, channel_info.xdim(0));
    EXPECT_EQ(2048, channel_info.ydim(0));
    EXPECT_EQ(3, channel_info.zdim(0));
    EXPECT_EQ(512, channel_info.xchunk_dim(0));
    EXPECT_EQ(256, channel_info.ychunk_dim(0));
    EXPECT_EQ(1, channel_info.zchunk_dim(0));
    EXPECT_EQ(2048, channel_info.slice_dim(0, mcp3d::ChannelAxis::X));
    EXPECT_EQ(2048, channel_info.slice_dim(0, mcp3d::ChannelAxis::Y));
    EXPECT_EQ(3, channel_info.slice_dim(0, mcp3d::ChannelAxis::Z));
    EXPECT_EQ(1, channel_info.slice_number(0, mcp3d::ChannelAxis::X));
    EXPECT_EQ(1, channel_info.slice_number(0, mcp3d::ChannelAxis::Y));
    EXPECT_EQ(1, channel_info.slice_number(0, mcp3d::ChannelAxis::Z));
    EXPECT_EQ(1, channel_info.TotalSliceNumber(0));
    EXPECT_EQ(96, channel_info.TotalNumberOfImageChunks(0));
    EXPECT_EQ(96, channel_info.ChunkNumberPerSlice(0));
    EXPECT_EQ(4, channel_info.AxisChunkNumberPerSlice(0, mcp3d::ChannelAxis::X));
    EXPECT_EQ(8, channel_info.AxisChunkNumberPerSlice(0, mcp3d::ChannelAxis::Y));
    EXPECT_EQ(3, channel_info.AxisChunkNumberPerSlice(0, mcp3d::ChannelAxis::Z));
    EXPECT_EQ(mcp3d::JoinPath(flat_dir, "flat_z000001_y000768_x000512.ome.tif"),
              channel_info.ImagePath(0, 1, 1000, 1000));
    // zyx slices (based on tiff/imaris data)
    // pyr_level_00: tif dim 291 * 6197 * 5183
    //               ome tif dim 320 * 8192 * 6144
    //               slice dim 16 * 2048 * 2048
    //               chunk dim 1 * 512 * 512
    string zyx_slice_dir(test_ome_tiff_format_dir("full_sagittal_section_300uM_cropped_challenging_ome_tiffs1"));
    mcp3d::MImage image(zyx_slice_dir);
    image.ReadImageInfo(0);
    EXPECT_EQ(5, image.channel_info(0).n_pyr_infos());
    EXPECT_EQ(0, image.channel_info(0).channel_number());
    EXPECT_EQ(1, image.image_info().n_channels(0));
    EXPECT_EQ(320, image.channel_info(0).zdim(0));
    EXPECT_EQ(8192, image.channel_info(0).ydim(0));
    EXPECT_EQ(6144, image.channel_info(0).xdim(0));
    EXPECT_EQ(16, image.channel_info(0).slice_dim(0, mcp3d::ChannelAxis::Z));
    EXPECT_EQ(2048, image.channel_info(0).slice_dim(0, mcp3d::ChannelAxis::Y));
    EXPECT_EQ(2048, image.channel_info(0).slice_dim(0, mcp3d::ChannelAxis::X));
    EXPECT_EQ(20, image.channel_info(0).slice_number(0, mcp3d::ChannelAxis::Z));
    EXPECT_EQ(4, image.channel_info(0).slice_number(0, mcp3d::ChannelAxis::Y));
    EXPECT_EQ(3, image.channel_info(0).slice_number(0, mcp3d::ChannelAxis::X));
    EXPECT_EQ(1, image.channel_info(0).zchunk_dim(0));
    EXPECT_EQ(512, image.channel_info(0).ychunk_dim(0));
    EXPECT_EQ(512, image.channel_info(0).xchunk_dim(0));

    image.SaveImageInfo();
}

TEST(MOmeTiffFormat, ReadChannelData)
{
    string imaris_dir(test_hdf5_format_dir()),
           tiff_dir(test_hdf5_format_dir("full_sagittal_section_300uM_cropped_challenging_tiff_planes")),
           ome_tiff_dir(test_ome_tiff_format_dir("full_sagittal_section_300uM_cropped_challenging_ome_tiffs1"));
    mcp3d::MImage image_imaris(imaris_dir), image_tiff(tiff_dir), image_ome_tiff(ome_tiff_dir);
    image_imaris.ReadImageInfo(0);
    image_tiff.ReadImageInfo(0);
    image_ome_tiff.ReadImageInfo(0);
    EXPECT_EQ(image_imaris.channel_info(0).zdim(0), image_tiff.channel_info(0).zdim(0));
    EXPECT_EQ(image_imaris.channel_info(0).ydim(0), image_tiff.channel_info(0).ydim(0));
    EXPECT_EQ(image_imaris.channel_info(0).xdim(0), image_tiff.channel_info(0).xdim(0));

    int64_t seed = std::chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);
    uniform_int_distribution<int> distribution_zoffset(0, 291),
            distribution_yoffset(0, 6196),
            distribution_xoffset(0, 5182);
    uniform_int_distribution<int> distribution_zextent(10, 50),
            distribution_yextent(512, 2048),
            distribution_xextent(512, 2048);
    uniform_int_distribution<int> distribution_zstride(1, 3),
            distribution_ystride(1, 5),
            distribution_xstride(1, 5);

    for (int i = 0; i < 5; ++i)
    {
        vector<int> offsets({distribution_zoffset(generator),
                             distribution_yoffset(generator),
                             distribution_xoffset(generator)});
        vector<int> extents({distribution_zextent(generator),
                             distribution_yextent(generator),
                             distribution_xextent(generator)});
        vector<int> strides({distribution_zstride(generator),
                             distribution_ystride(generator),
                             distribution_xstride(generator)});
        mcp3d::MImageBlock block(offsets, extents, strides);
        INIT_TIMER(0)
        INIT_TIMER(1)
        INIT_TIMER(2)
        image_imaris.SelectView(block, 0);
        TIC(0)
        image_imaris.ReadData();
        TOC(0)
        image_tiff.SelectView(block, 0);
        TIC(1)
        image_tiff.ReadData();
        TOC(1)
        image_ome_tiff.SelectView(block, 0);
        TIC(2)
        image_ome_tiff.ReadData();
        TOC(2)
        EXPECT_TRUE(image_ome_tiff.HasEqualData(image_imaris));
        EXPECT_TRUE(image_ome_tiff.HasEqualData(image_tiff));
        REPORT_TIME_TO_COMPLETION("time to read from imaris file", 0)
        REPORT_TIME_TO_COMPLETION("time to read from tiff sequence", 1)
        REPORT_TIME_TO_COMPLETION("time to read from ome tiff", 2)
    }

    //mcp3d::RemovePath(image_ome_tiff.channel_info(0).channel_info_path());
}