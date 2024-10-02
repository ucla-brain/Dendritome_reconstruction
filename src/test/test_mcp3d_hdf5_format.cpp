//
// Created by muyezhu on 3/12/19.
//
#include <tuple>
#include <vector>
#include <random>
#include <chrono>
#include <hdf5.h>
#include <gtest/gtest.h>
#include "common/mcp3d_common.hpp"
#include "image_interface/mcp3d_voxel_types.hpp"
#include "image_interface/mcp3d_hdf5_utils.hpp"
#include "image_interface/mcp3d_imaris_util.hpp"
#include "image/mcp3d_channel_info.hpp"
#include "image/mcp3d_image_io.hpp"

using namespace std;

string test_hdf5_format_dir(const string& subdir = string())
{
    string d = mcp3d::JoinPath({mcp3d::test_data_dir(),
                                "image_formats", "hdf5_format", subdir});
    MCP3D_ASSERT(mcp3d::IsDir(d))
    return d;
}

TEST(MImarisFormat, ImarisResolutions)
{
    string imaris_path(mcp3d::JoinPath({test_hdf5_format_dir(),
                                        "full_sagittal_section_300uM_cropped_challenging_FusionStitcher.ims"}));
    EXPECT_TRUE(mcp3d::IsFile(imaris_path));
    vector<string> results = mcp3d::ImarisResolutionNames(imaris_path);
    vector<string> groups({"ResolutionLevel 0", "ResolutionLevel 1",
                           "ResolutionLevel 2", "ResolutionLevel 3",
                           "ResolutionLevel 4", "ResolutionLevel 5",
                           "ResolutionLevel 6"});
    EXPECT_EQ(groups, results);
}

TEST(MImarisFormat, ImarisResolutionChannels)
{
    string imaris_path(mcp3d::JoinPath({test_hdf5_format_dir(),
                                        "2019-01-18_15.51.15_Protocol_F00.ims"}));
    EXPECT_TRUE(mcp3d::IsFile(imaris_path));
    hid_t imaris_id = mcp3d::Hdf5Handle(imaris_path);
    vector<string> results = mcp3d::ImarisResolutionChannelNames(imaris_id, 0);
    vector<string> channels({"Channel 0", "Channel 1", "Channel 2", "Channel 3"});
    EXPECT_EQ(channels, results);
    results = mcp3d::ImarisResolutionChannelNames(imaris_id, 1);
    EXPECT_EQ(channels, results);
    H5Fclose(imaris_id);
}

TEST(MImarisFormat, ImarisImageXyzSizes)
{
    string imaris_path(mcp3d::JoinPath({test_hdf5_format_dir(),
                                        "2019-01-18_15.51.15_Protocol_F00.ims"}));
    EXPECT_TRUE(mcp3d::IsFile(imaris_path));
    hid_t imaris_id = mcp3d::Hdf5Handle(imaris_path);
    vector<int> xyz_dims = mcp3d::ImarisChannelImageXyzSizes(imaris_id, 0, 0,
                                                             0);
    EXPECT_EQ(2, xyz_dims[0]);
    EXPECT_EQ(2048, xyz_dims[1]);
    EXPECT_EQ(2048, xyz_dims[2]);
    H5Fclose(imaris_id);

    imaris_path = mcp3d::JoinPath({test_hdf5_format_dir(),
                                   "full_sagittal_section_300uM_cropped_challenging_FusionStitcher.ims"});
    EXPECT_TRUE(mcp3d::IsFile(imaris_path));
    imaris_id = mcp3d::Hdf5Handle(imaris_path);
    xyz_dims = mcp3d::ImarisChannelImageXyzSizes(imaris_id, 0, 0, 0);
    EXPECT_EQ(292, xyz_dims[0]);
    EXPECT_EQ(6197, xyz_dims[1]);
    EXPECT_EQ(5183, xyz_dims[2]);
    H5Fclose(imaris_id);
}

TEST(MImarisFormat, Hdf5DatasetChunkDimensions)
{
    string imaris_path = mcp3d::JoinPath({test_hdf5_format_dir(),
                                          "full_sagittal_section_300uM_cropped_challenging_FusionStitcher.ims"});
    mcp3d::ImarisResolutionLevelInfo info(imaris_path, 0, 0);
    EXPECT_EQ(vector<int>({8, 256, 256}), info.chunk_xyz_dims());
}

TEST(MImarisFormat, Hdf5DatasetVoxelType)
{
    string imaris_path(mcp3d::JoinPath({test_hdf5_format_dir(), "2019-01-18_15.51.15_Protocol_F00.ims"}));
    EXPECT_TRUE(mcp3d::IsFile(imaris_path));
    hid_t imaris_id = mcp3d::Hdf5Handle(imaris_path);
    EXPECT_EQ(mcp3d::VoxelType::M16U,
              mcp3d::ImarisResolutionVoxelType(imaris_id, 0, 0));
}

void MakeMockHdf5File()
{

}

TEST(MImarisFormat, SetHdf5DatasetZlibDeflate)
{
    // zlib filter is supported by hdf5 library
    EXPECT_TRUE(mcp3d::Hdf5ZlibFilterAvailable());

}

TEST(MImarisFormat, ReadChannelInfo)
{
    mcp3d::MImageIO image_io {};
    mcp3d::MImageInfo image_info(test_hdf5_format_dir(), {""});
    mcp3d::MChannelInfo channel_info = image_io.ReadChannelInfo(image_info, 0);
    // ResolutionLevel 6 will be discarded due to insufficient z dimension
    EXPECT_EQ(6, channel_info.n_pyr_infos());
    EXPECT_EQ(mcp3d::FileFormat::IMARIS, channel_info.file_format(0));
    EXPECT_EQ(mcp3d::VoxelType::M16U, channel_info.voxel_type(0));
    EXPECT_EQ(test_hdf5_format_dir(), channel_info.channel_root_dir());
    EXPECT_EQ(vector<int>({292, 6197, 5183}), channel_info.xyz_dims(0));
    EXPECT_EQ(test_hdf5_format_dir(), channel_info.channel_root_dir());
    EXPECT_EQ(vector<int>({1, 2, 4, 8, 16, 32}),
              channel_info.pyr_xy_ratios());
    EXPECT_EQ(vector<int>({1, 1, 2, 4, 8, 16}),
              channel_info.pyr_z_ratios());
    // make pyr_level_01 and expect error (image root is composite,
    // should have no more composite levels)
    string target(mcp3d::JoinPath({test_hdf5_format_dir(),
                                   "full_sagittal_section_300uM_cropped_challenging_tiff_planes/pyr_level_01"}));
    string link(mcp3d::JoinPath({test_hdf5_format_dir(), "pyr_level_01"}));
    string cmd("ln -s ");
    cmd.append(target);
    cmd.append(" ");
    cmd.append(link);
    mcp3d::SysCmdResult(cmd);
    EXPECT_TRUE(mcp3d::IsDir(link));
    EXPECT_THROW(channel_info = image_io.ReadChannelInfo(image_info, 0),
                 mcp3d::MCPAssertionError);
    cmd = "rm ";
    cmd.append(link);
    mcp3d::SysCmdResult(cmd);
    EXPECT_FALSE(mcp3d::IsDir(link));
}

TEST(MImarisFormat, ReadDataImaris)
{
    mcp3d::MImage image_imaris(test_hdf5_format_dir()),
                        image_tiff(mcp3d::JoinPath({test_hdf5_format_dir(),
                                                    "full_sagittal_section_300uM_cropped_challenging_tiff_planes"}));
    // 292 * 6197 * 5183
    image_imaris.ReadImageInfo(0);
    // the tiff hierarchies are read directly from .ims file and written
    // as tif sequences. xyz dimensions are ImageSizeX, ImageSizeY, ImageSizeZ
    image_tiff.ReadImageInfo(0);
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
        image_imaris.SelectView(block, 0);
        TIC(0)
        image_imaris.ReadData();
        TOC(0)
        image_tiff.SelectView(block, 0);
        TIC(1)
        image_tiff.ReadData();
        TOC(1)
        EXPECT_TRUE(image_imaris.HasEqualData(image_tiff));
        REPORT_TIME_TO_COMPLETION("time to read from imaris file", 0)
        REPORT_TIME_TO_COMPLETION("time to read from tiff sequence", 1)
    }

}

