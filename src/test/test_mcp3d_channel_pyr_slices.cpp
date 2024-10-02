//
// Created by muyezhu on 9/13/20.
//
#include <fstream>
#include <tuple>
#include <unordered_map>
#include <regex>
#include <boost/algorithm/string/predicate.hpp>
#include "common/mcp3d_common.hpp"
#include "image_interface/mcp3d_file_formats.hpp"
#include "image_layout/mcp3d_channel_pyr_slices.hpp"
#include "test_mcp3d_channel_pyr_slices.hpp"

using namespace std;

void MChannelPyrSlicesTest::Configure(const string &dir_name, const vector<int> &slice_numbers, const vector<int> &slice_dims, bool dir_name_is_full_path)
{
    // remove previously made directory
    CleanUp();
    channel_pyr_dir_ = dir_name_is_full_path ? dir_name : mcp3d::JoinPath(test_channel_pyr_slices_dir(), dir_name);
    mcp3d::RemovePath(channel_pyr_dir_);
    ASSERT_FALSE(mcp3d::IsDir(channel_pyr_dir_));
    mcp3d::MakeDirectories(channel_pyr_dir_);
    ASSERT_TRUE(mcp3d::IsDir(channel_pyr_dir_));
    slice_numbers_ = slice_numbers;
    slice_dims_ = slice_dims;
    for (int i = 0; i < 3; ++i)
        ASSERT_EQ(slice_numbers_[i] == 1, slice_dims_[i] == mcp3d::MAX_AXIS_SLICE_DIM);

    string slice_dir;
    for (int zid = 0; zid < slice_numbers[0]; ++zid)
        for (int yid = 0; yid < slice_numbers[1]; ++yid)
            for (int xid = 0; xid < slice_numbers[2]; ++xid)
            {
                GetSliceDir(zid, yid, xid, slice_dir);
                mcp3d::MakeDirectories(slice_dir);
                ASSERT_TRUE(mcp3d::IsDir(slice_dir));
            }
    pyr_slices_.set_channel_pyr_dir(channel_pyr_dir_);
}

void MChannelPyrSlicesTest::AddDirs(const vector<string> &dir_paths)
{
    for (const auto& dir_path: dir_paths)
    {
        ASSERT_TRUE(boost::algorithm::starts_with(dir_path, channel_pyr_dir_));
        mcp3d::MakeDirectories(dir_path);
    }
    pyr_slices_.Refresh();
}

void MChannelPyrSlicesTest::AddFiles(const vector<string> &file_names, int slice_zid, int slice_yid, int slice_xid)
{
    string slice_dir;
    GetSliceDir(slice_zid, slice_yid, slice_xid, slice_dir);
    for (const auto& file_name: file_names)
    {
        ofstream of{mcp3d::JoinPath(slice_dir, file_name)};
        of.close();
    }
    pyr_slices_.Refresh();
}

void MChannelPyrSlicesTest::RemoveDir(int slice_zid, int slice_yid, int slice_xid)
{
    vector<string> axis_dirs;
    vector<mcp3d::ChannelAxis> axes({mcp3d::ChannelAxis::Z, mcp3d::ChannelAxis::Y, mcp3d::ChannelAxis::X});
    vector<int> slice_ids({slice_zid, slice_yid, slice_xid});
    for (int i = 0; i < 3; ++i)
    {
        if (slice_ids[i] < 0)
            break;
        ASSERT_TRUE(slice_ids[i] < slice_numbers_[i]);
        axis_dirs.push_back(mcp3d::MChannelPyrSlices::AxisSliceName(axes[i], slice_ids[i], slice_dims_[i]));
    }
    string dir_name = mcp3d::JoinPath(axis_dirs);
    mcp3d::RemovePath(mcp3d::JoinPath(channel_pyr_dir_, dir_name));
    pyr_slices_.Refresh();
}

void MChannelPyrSlicesTest::RemoveFiles(const vector<string> &file_names, int slice_zid, int slice_yid, int slice_xid)
{
    string slice_dir;
    GetSliceDir(slice_zid, slice_yid, slice_xid, slice_dir);
    for (const auto& file_name: file_names)
        mcp3d::RemovePath(mcp3d::JoinPath(slice_dir, file_name));
    pyr_slices_.Refresh();
}

void MChannelPyrSlicesTest::RefreshPyrSlices()
{
    pyr_slices_.Refresh();
}

void MChannelPyrSlicesTest::GetSliceDir(int slice_zid, int slice_yid, int slice_xid, string& slice_dir) const
{
    /// when ASSERT_* is used outside of TEST macro, the function needs to be of void return type. otherwise doesn't compile
    ASSERT_TRUE(mcp3d::IsDir(channel_pyr_dir_));
    ASSERT_TRUE(slice_zid >= 0 && slice_zid < slice_numbers_[0]);
    ASSERT_TRUE(slice_yid >= 0 && slice_yid < slice_numbers_[1]);
    ASSERT_TRUE(slice_xid >= 0 && slice_xid < slice_numbers_[2]);
    vector<mcp3d::ChannelAxis> axes({mcp3d::ChannelAxis::Z, mcp3d::ChannelAxis::Y, mcp3d::ChannelAxis::X});
    slice_dir = mcp3d::JoinPath(channel_pyr_dir_, mcp3d::MChannelPyrSlices::AxisSliceName(axes[0], slice_zid, slice_dims_[0]),
                                mcp3d::MChannelPyrSlices::AxisSliceName(axes[1], slice_yid, slice_dims_[1]),
                                mcp3d::MChannelPyrSlices::AxisSliceName(axes[2], slice_xid, slice_dims_[2]));
}

TEST(MChannelPyrSlices, AxisSliceNamePattern)
{
    smatch sm;
    string input = "x0000_123456";
    EXPECT_TRUE(regex_match(input, sm, mcp3d::MChannelPyrSlices::AxisSliceNamePattern()));
    EXPECT_EQ("x", sm.str(1));
    EXPECT_EQ(0, stoi(sm.str(2)));
    EXPECT_EQ(123456, stoi(sm.str(3)));
    input = "y9212_000621";
    EXPECT_TRUE(regex_match(input, sm, mcp3d::MChannelPyrSlices::AxisSliceNamePattern()));
    EXPECT_EQ("y", sm.str(1));
    EXPECT_EQ(9212, stoi(sm.str(2)));
    EXPECT_EQ(621, stoi(sm.str(3)));
    input = string("z0023_000008");
    EXPECT_TRUE(regex_match(input, sm, mcp3d::MChannelPyrSlices::AxisSliceNamePattern()));
    EXPECT_EQ("z", sm.str(1));
    EXPECT_EQ(23, stoi(sm.str(2)));
    EXPECT_EQ(8, stoi(sm.str(3)));

    input = "";
    EXPECT_FALSE(regex_match(input, mcp3d::MChannelPyrSlices::AxisSliceNamePattern()));
    input = "X0000_123456";
    EXPECT_FALSE(regex_match(input, mcp3d::MChannelPyrSlices::AxisSliceNamePattern()));
    input = "x000_123456";
    EXPECT_FALSE(regex_match(input, mcp3d::MChannelPyrSlices::AxisSliceNamePattern()));
    input = "x0000_12345";
    EXPECT_FALSE(regex_match(input, mcp3d::MChannelPyrSlices::AxisSliceNamePattern()));
    input = "x0000_a12345";
    EXPECT_FALSE(regex_match(input, mcp3d::MChannelPyrSlices::AxisSliceNamePattern()));
    input = "yx0000_12345";
    EXPECT_FALSE(regex_match(input, mcp3d::MChannelPyrSlices::AxisSliceNamePattern()));
    input = "x0000_123456y";
    EXPECT_FALSE(regex_match(input, mcp3d::MChannelPyrSlices::AxisSliceNamePattern()));
}

TEST(MChannelPyrSlices, AxisSliceName)
{
    EXPECT_EQ("", mcp3d::MChannelPyrSlices::AxisSliceName(mcp3d::ChannelAxis::X, 0, mcp3d::MAX_AXIS_SLICE_DIM));
    mcp3d::ChannelAxis zaxis = mcp3d::ChannelAxis::Z;
    EXPECT_EQ("z0003_000256", mcp3d::MChannelPyrSlices::AxisSliceName(zaxis, 3, 256));
    EXPECT_THROW(mcp3d::MChannelPyrSlices::AxisSliceName(zaxis, -3, 256), mcp3d::MCPAssertionError);
    EXPECT_THROW(mcp3d::MChannelPyrSlices::AxisSliceName(zaxis, 30000, 256), mcp3d::MCPAssertionError);
    EXPECT_THROW(mcp3d::MChannelPyrSlices::AxisSliceName(zaxis, 3, -256), mcp3d::MCPAssertionError);
    EXPECT_THROW(mcp3d::MChannelPyrSlices::AxisSliceName(zaxis, 3, 2560000), mcp3d::MCPAssertionError);
}

TEST(MChannelPyrSlices, Constructor)
{
    mcp3d::MChannelPyrSlices empty;
    EXPECT_TRUE(empty.empty());
    EXPECT_TRUE(empty.channel_pyr_dir().empty());
    EXPECT_EQ(mcp3d::FileFormat::UNKNOWN, empty.file_format());
    for (const auto& axis: vector<mcp3d::ChannelAxis>({mcp3d::ChannelAxis::Z, mcp3d::ChannelAxis::Y, mcp3d::ChannelAxis::X}))
    {
        EXPECT_EQ(0, empty.slice_number(axis));
        EXPECT_EQ(-1, empty.slice_dim(axis));
    }
    string bad_dir = "bad_dir!";
    EXPECT_THROW(mcp3d::MChannelPyrSlices{bad_dir}, mcp3d::MCPOSError);
}

TEST_F(MChannelPyrSlicesTest, ReadChannelPyrSlices)
{
    // (dir name, slice numbers, slice dims)
    vector<tuple<string, vector<int>, vector<int>>> configs =
    {
        tuple<string, vector<int>, vector<int>>("flat", vector<int>(3, 1), vector<int>(3, mcp3d::MAX_AXIS_SLICE_DIM)),
        tuple<string, vector<int>, vector<int>>("full", vector<int>({4, 5, 6}), vector<int>({1024, 2048, 256})),
        tuple<string, vector<int>, vector<int>>("half_full0", vector<int>({4, 1, 6}), vector<int>({1024, mcp3d::MAX_AXIS_SLICE_DIM, 256})),
        tuple<string, vector<int>, vector<int>>("half_full1", vector<int>({4, 1, 1}), vector<int>({1024, mcp3d::MAX_AXIS_SLICE_DIM, mcp3d::MAX_AXIS_SLICE_DIM})),
        tuple<string, vector<int>, vector<int>>("half_full2", vector<int>({1, 1, 6}), vector<int>({mcp3d::MAX_AXIS_SLICE_DIM, mcp3d::MAX_AXIS_SLICE_DIM, 256}))
    };

    for (const auto& config: configs)
    {
        Configure(get<0>(config), get<1>(config), get<2>(config));
        EXPECT_FALSE(pyr_slices_.empty());
        EXPECT_EQ(get<1>(config)[0], pyr_slices_.slice_number(mcp3d::ChannelAxis::Z));
        EXPECT_EQ(get<1>(config)[1], pyr_slices_.slice_number(mcp3d::ChannelAxis::Y));
        EXPECT_EQ(get<1>(config)[2], pyr_slices_.slice_number(mcp3d::ChannelAxis::X));
        EXPECT_EQ(get<2>(config)[0], pyr_slices_.slice_dim(mcp3d::ChannelAxis::Z));
        EXPECT_EQ(get<2>(config)[1], pyr_slices_.slice_dim(mcp3d::ChannelAxis::Y));
        EXPECT_EQ(get<2>(config)[2], pyr_slices_.slice_dim(mcp3d::ChannelAxis::X));
        EXPECT_EQ(mcp3d::JoinPath(test_channel_pyr_slices_dir(), get<0>(config)), pyr_slices_.channel_pyr_dir());
        EXPECT_EQ(mcp3d::FileFormat::UNKNOWN, pyr_slices_.file_format());
    }
    // full config. incomplete slice directories
    Configure(get<0>(configs[1]), get<1>(configs[1]), get<2>(configs[1]));
    try
    {
        AddDirs({mcp3d::JoinPath(channel_pyr_dir_, mcp3d::MChannelPyrSlices::AxisSliceName(mcp3d::ChannelAxis::Z, 0, slice_dims_[0]),
                                 mcp3d::MChannelPyrSlices::AxisSliceName(mcp3d::ChannelAxis::Y, 0, slice_dims_[1]),
                                 mcp3d::MChannelPyrSlices::AxisSliceName(mcp3d::ChannelAxis::X, slice_numbers_[2], slice_dims_[2]))});
    }
    catch (const mcp3d::MCPRuntimeError& e)
    {
        EXPECT_TRUE(boost::algorithm::ends_with(string(e.what()), "expected slice directory " +
                                                                  mcp3d::JoinPath(channel_pyr_dir_, "z0000_001024/y0001_002048/x0006_000256") + " not found"));
    }
    catch (...)
    {
        ADD_FAILURE();
    }
    // half full1 config. single slice along an axis has dimension not equal to MAX_AXIS_SLICE_DIM
    Configure(get<0>(configs[3]), get<1>(configs[3]), get<2>(configs[3]));
    try
    {
        AddDirs({mcp3d::JoinPath(channel_pyr_dir_, "z0000_001024", "y0000_000256")});
    }
    catch (const mcp3d::MCPRuntimeError& e)
    {
        EXPECT_TRUE(boost::algorithm::ends_with(string(e.what()), "axis with a single slice should always have "
                                                                  "mcp3d::MAX_AXIS_SLICE_DIM as slice dimension, and vice versa"));
    }
    catch (...)
    {
        ADD_FAILURE();
    }
    // half full1 config. slices from mixed axes under same level
    Configure(get<0>(configs[3]), get<1>(configs[3]), get<2>(configs[3]));
    try
    {
        AddDirs({mcp3d::JoinPath(channel_pyr_dir_, "y0000_000256")});
    }
    catch (const mcp3d::MCPRuntimeError& e)
    {
        EXPECT_TRUE(boost::algorithm::ends_with(string(e.what()), "folder slices under " + channel_pyr_dir_ + " do not belong to the same axis"));
    }
    catch (...)
    {
        ADD_FAILURE();
    }
    // half full1 config. mixed slice dim values
    Configure(get<0>(configs[3]), get<1>(configs[3]), get<2>(configs[3]));
    EXPECT_THROW(AddDirs({mcp3d::JoinPath(channel_pyr_dir_, "z0000_001025")}), mcp3d::MCPAssertionError);
    // half full1 config. zyx axis ordering
    Configure(get<0>(configs[3]), get<1>(configs[3]), get<2>(configs[3]));
    try
    {
        AddDirs({mcp3d::JoinPath(channel_pyr_dir_, "z0000_001024", "z0000_001024")});
    }
    catch (const mcp3d::MCPRuntimeError& e)
    {
        EXPECT_TRUE(boost::algorithm::ends_with(string(e.what()), "folder slices under " + mcp3d::JoinPath(channel_pyr_dir_, "z0000_001024") +
                                                                  " do not follow the expected zyx axes ordering"));
    }
    catch (...)
    {
        ADD_FAILURE();
    }
}

TEST_F(MChannelPyrSlicesTest, ReadChannelPyrFormat)
{
    // (dir name, slice numbers, slice dims)
    vector<tuple<string, vector<int>, vector<int>>> configs =
    {
        // imaris, hdf5, tif, ome.tif
        tuple<string, vector<int>, vector<int>>("flat", vector<int>(3, 1), vector<int>(3, mcp3d::MAX_AXIS_SLICE_DIM)),
        // hdf5, ome.tif
        tuple<string, vector<int>, vector<int>>("full", vector<int>({4, 5, 6}), vector<int>({1024, 2048, 256})),
        // hdf5, tif, ome.tif
        tuple<string, vector<int>, vector<int>>("half_full", vector<int>({4, 1, 1}), vector<int>({1024, mcp3d::MAX_AXIS_SLICE_DIM, mcp3d::MAX_AXIS_SLICE_DIM}))
    };

    // all formats compatible
    Configure(get<0>(configs[0]), get<1>(configs[0]), get<2>(configs[0]));
    // (1) add 1 stitched imaris file
    AddFiles(vector<string>({"1.FusionStitcher.ims", "1.hdf5", "1.tif", "1.ome.tif", "1.ims", "1.fake"}), 0, 0, 0);
    EXPECT_EQ(mcp3d::FileFormat::IMARIS, pyr_slices_.file_format());
    // (2) add 1 stitched imaris file
    EXPECT_THROW(AddFiles(vector<string>(1, "2.FusionStitcher.ims"), 0, 0, 0), mcp3d::MCPRuntimeError);
    // (3) remove 1 stitched imaris file
    RemoveFiles(vector<string>(1, "2.FusionStitcher.ims"), 0, 0, 0);
    EXPECT_EQ(mcp3d::FileFormat::IMARIS, pyr_slices_.file_format());
    // (4) remove 1 stitched imaris file
    RemoveFiles(vector<string>(1, "1.FusionStitcher.ims"), 0, 0, 0);
    EXPECT_EQ(mcp3d::FileFormat::HDF5, pyr_slices_.file_format());
    // (5) remove 1 hdf5 file
    RemoveFiles(vector<string>(1, "1.hdf5"), 0, 0, 0);
    EXPECT_EQ(mcp3d::FileFormat::TIFF, pyr_slices_.file_format());
    // (6) remove 1 tiff file
    RemoveFiles(vector<string>(1, "1.tif"), 0, 0, 0);
    EXPECT_EQ(mcp3d::FileFormat::OMETIFF, pyr_slices_.file_format());
    // (7) remove 1 ome.tiff file
    RemoveFiles(vector<string>(1, "1.ome.tif"), 0, 0, 0);
    EXPECT_EQ(mcp3d::FileFormat::UNKNOWN, pyr_slices_.file_format());

    // hdf5 and ome tiff compatible
    Configure(get<0>(configs[1]), get<1>(configs[1]), get<2>(configs[1]));
    // add 1 stitched imaris file and 1 tiff file in slice (0, 0, 0) (formats not compatible with directory)
    // add 1 ome.tif file in (0, 0, 1), 1 hdf file in (0, 0, 2)
    AddFiles(vector<string>({"1.FusionStitcher.ims", "1.tif", "1.fake"}), 0, 0, 0);
    EXPECT_EQ(mcp3d::FileFormat::UNKNOWN, pyr_slices_.file_format());
    AddFiles(vector<string>(1, "1.ome.tif"), 0, 0, 1);
    EXPECT_EQ(mcp3d::FileFormat::OMETIFF, pyr_slices_.file_format());
    AddFiles(vector<string>(1, "1.hdf5"), 0, 0, 2);
    EXPECT_EQ(mcp3d::FileFormat::HDF5, pyr_slices_.file_format());

    // hdf5, tif and ome tiff compatible
    Configure(get<0>(configs[2]), get<1>(configs[2]), get<2>(configs[2]));
    AddFiles(vector<string>({"1.FusionStitcher.ims", "1.tif", "1.fake"}), 0, 0, 0);
    EXPECT_EQ(mcp3d::FileFormat::TIFF, pyr_slices_.file_format());
    AddFiles(vector<string>({"1.ome.tif", "1.hdf5", "1.tif", "1.fake"}), 0, 0, 0);
    EXPECT_EQ(mcp3d::FileFormat::HDF5, pyr_slices_.file_format());
}


TEST_F(MChannelPyrSlicesTest, Clear)
{
    tuple<string, vector<int>, vector<int>> config("half_full", vector<int>({4, 2, 1}), vector<int>({1024, 256, mcp3d::MAX_AXIS_SLICE_DIM}));
    Configure(get<0>(config), get<1>(config), get<2>(config));
    AddFiles(vector<string>({"1.hdf5"}), 0, 0, 0);
    EXPECT_FALSE(pyr_slices_.empty());
    EXPECT_EQ(get<1>(config)[0], pyr_slices_.slice_number(mcp3d::ChannelAxis::Z));
    EXPECT_EQ(get<1>(config)[1], pyr_slices_.slice_number(mcp3d::ChannelAxis::Y));
    EXPECT_EQ(get<1>(config)[2], pyr_slices_.slice_number(mcp3d::ChannelAxis::X));
    EXPECT_EQ(get<2>(config)[0], pyr_slices_.slice_dim(mcp3d::ChannelAxis::Z));
    EXPECT_EQ(get<2>(config)[1], pyr_slices_.slice_dim(mcp3d::ChannelAxis::Y));
    EXPECT_EQ(get<2>(config)[2], pyr_slices_.slice_dim(mcp3d::ChannelAxis::X));
    EXPECT_EQ(channel_pyr_dir_, pyr_slices_.channel_pyr_dir());
    EXPECT_EQ(mcp3d::FileFormat::HDF5, pyr_slices_.file_format());

    pyr_slices_.Clear();
    EXPECT_TRUE(pyr_slices_.empty());
    for (const auto& axis: vector<mcp3d::ChannelAxis>({mcp3d::ChannelAxis::Z, mcp3d::ChannelAxis::Y, mcp3d::ChannelAxis::X}))
    {
        EXPECT_EQ(0, pyr_slices_.slice_number(axis));
        EXPECT_EQ(-1, pyr_slices_.slice_dim(axis));
    }
    EXPECT_TRUE(pyr_slices_.channel_pyr_dir().empty());
    EXPECT_EQ(mcp3d::FileFormat::UNKNOWN, pyr_slices_.file_format());
}

TEST_F(MChannelPyrSlicesTest, Refresh)
{
    tuple<string, vector<int>, vector<int>> config("half_full", vector<int>({4, 2, 1}), vector<int>({1024, 256, mcp3d::MAX_AXIS_SLICE_DIM}));
    Configure(get<0>(config), get<1>(config), get<2>(config));
    AddFiles(vector<string>({"1.hdf5"}), 0, 0, 0);
    EXPECT_FALSE(pyr_slices_.empty());
    EXPECT_EQ(get<1>(config)[0], pyr_slices_.slice_number(mcp3d::ChannelAxis::Z));
    EXPECT_EQ(get<1>(config)[1], pyr_slices_.slice_number(mcp3d::ChannelAxis::Y));
    EXPECT_EQ(get<1>(config)[2], pyr_slices_.slice_number(mcp3d::ChannelAxis::X));
    EXPECT_EQ(get<2>(config)[0], pyr_slices_.slice_dim(mcp3d::ChannelAxis::Z));
    EXPECT_EQ(get<2>(config)[1], pyr_slices_.slice_dim(mcp3d::ChannelAxis::Y));
    EXPECT_EQ(get<2>(config)[2], pyr_slices_.slice_dim(mcp3d::ChannelAxis::X));
    EXPECT_EQ(channel_pyr_dir_, pyr_slices_.channel_pyr_dir());
    EXPECT_EQ(mcp3d::FileFormat::HDF5, pyr_slices_.file_format());

    EXPECT_TRUE(mcp3d::RemovePath(channel_pyr_dir_));

}

TEST_F(MChannelPyrSlicesTest, AxisIsFlat)
{
    mcp3d::MChannelPyrSlices empty;
    EXPECT_FALSE(empty.IsFlat());
    EXPECT_FALSE(empty.AxisIsFlat(mcp3d::ChannelAxis::Z));
    EXPECT_FALSE(empty.AxisIsFlat(mcp3d::ChannelAxis::Y));
    EXPECT_FALSE(empty.AxisIsFlat(mcp3d::ChannelAxis::X));

    // (dir name, slice numbers, slice dims)
    vector<tuple<string, vector<int>, vector<int>>> configs =
    {
        tuple<string, vector<int>, vector<int>>("flat", vector<int>(3, 1), vector<int>(3, mcp3d::MAX_AXIS_SLICE_DIM)),
        tuple<string, vector<int>, vector<int>>("half_full", vector<int>({4, 2, 1}), vector<int>({1024, 256, mcp3d::MAX_AXIS_SLICE_DIM}))
    };
    for (const auto& config: configs)
    {
        Configure(get<0>(config), get<1>(config), get<2>(config));
        EXPECT_EQ(get<1>(config)[0] == 1, pyr_slices_.AxisIsFlat(mcp3d::ChannelAxis::Z));
        EXPECT_EQ(get<1>(config)[1] == 1, pyr_slices_.AxisIsFlat(mcp3d::ChannelAxis::Y));
        EXPECT_EQ(get<1>(config)[2] == 1, pyr_slices_.AxisIsFlat(mcp3d::ChannelAxis::X));
        EXPECT_EQ(get<0>(config) == "flat", pyr_slices_.IsFlat());
    }
}

TEST_F(MChannelPyrSlicesTest, SliceName)
{
    mcp3d::MChannelPyrSlices empty;
    EXPECT_EQ("", empty.SliceNameFromSliceIDs(10, 11, 12));
    EXPECT_EQ("", empty.SliceNameFromCoordinates(-1, 0, 999));
    EXPECT_EQ(vector<string>{}, empty.SliceNames());

    tuple<string, vector<int>, vector<int>> config("half_full", vector<int>({4, 2, 1}), vector<int>({1024, 256, mcp3d::MAX_AXIS_SLICE_DIM}));
    Configure(get<0>(config), get<1>(config), get<2>(config));
    EXPECT_EQ("z0003_001024/y0001_000256", pyr_slices_.SliceNameFromSliceIDs(3, 1, 0));
    EXPECT_THROW(pyr_slices_.SliceNameFromSliceIDs(-3, 1, 0), mcp3d::MCPAssertionError);
    EXPECT_THROW(pyr_slices_.SliceNameFromSliceIDs(4, 1, 0), mcp3d::MCPAssertionError);
    EXPECT_EQ("z0000_001024/y0001_000256", pyr_slices_.SliceNameFromCoordinates(1000, 500, 99999));
    EXPECT_THROW(pyr_slices_.SliceNameFromCoordinates(1000, 500, -99999), mcp3d::MCPAssertionError);
    EXPECT_THROW(pyr_slices_.SliceNameFromCoordinates(1000, 512, 99999), mcp3d::MCPAssertionError);
}

TEST_F(MChannelPyrSlicesTest, CopyAndEquality)
{
    vector<tuple<string, vector<int>, vector<int>>> configs =
    {
        tuple<string, vector<int>, vector<int>>("half_full0", vector<int>({4, 3, 1}), vector<int>({1024, 256, mcp3d::MAX_AXIS_SLICE_DIM})),
        tuple<string, vector<int>, vector<int>>("half_full1", vector<int>({4, 2, 1}), vector<int>({1024, 256, mcp3d::MAX_AXIS_SLICE_DIM}))
    };

    Configure(get<0>(configs[0]), get<1>(configs[0]), get<2>(configs[0]));
    AddFiles(vector<string>({"1.FusionStitcher.ims"}), 0, 0, 0);

    mcp3d::MChannelPyrSlices copy0{pyr_slices_};
    EXPECT_EQ(pyr_slices_.slice_number(mcp3d::ChannelAxis::Z), copy0.slice_number(mcp3d::ChannelAxis::Z));
    EXPECT_EQ(pyr_slices_.slice_number(mcp3d::ChannelAxis::Y), copy0.slice_number(mcp3d::ChannelAxis::Y));
    EXPECT_EQ(pyr_slices_.slice_number(mcp3d::ChannelAxis::X), copy0.slice_number(mcp3d::ChannelAxis::X));
    EXPECT_EQ(pyr_slices_.slice_dim(mcp3d::ChannelAxis::Z), copy0.slice_dim(mcp3d::ChannelAxis::Z));
    EXPECT_EQ(pyr_slices_.slice_dim(mcp3d::ChannelAxis::Y), copy0.slice_dim(mcp3d::ChannelAxis::Y));
    EXPECT_EQ(pyr_slices_.slice_dim(mcp3d::ChannelAxis::X), copy0.slice_dim(mcp3d::ChannelAxis::X));
    EXPECT_EQ(pyr_slices_.channel_pyr_dir(), copy0.channel_pyr_dir());
    EXPECT_EQ(pyr_slices_.file_format(), copy0.file_format());
    EXPECT_TRUE(copy0 == pyr_slices_);
    EXPECT_FALSE(copy0 != pyr_slices_);

    Configure(get<0>(configs[1]), get<1>(configs[1]), get<2>(configs[1]));
    AddFiles(vector<string>({"1.FusionStitcher.ims"}), 0, 0, 0);
    mcp3d::MChannelPyrSlices copy1 = pyr_slices_;
    EXPECT_TRUE(copy1 == pyr_slices_);
    EXPECT_FALSE(copy1 != pyr_slices_);
    EXPECT_FALSE(copy0 == pyr_slices_);
    EXPECT_TRUE(copy0 != copy1);
}

