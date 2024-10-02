//
// Created by muyezhu on 3/24/18.
//
#include <cstring>
#include <random>
#include <algorithm>
#include <unordered_set>
#include <gtest/gtest.h>
#include "image_interface/mcp3d_file_formats.hpp"
#include "image_interface/mcp3d_hdf5_utils.hpp"
#include "image_interface/mcp3d_imaris_util.hpp"
#include "image_layout/mcp3d_image_layout_constants.hpp"
#include "image_layout/mcp3d_channel_layout.hpp"
#include "test_mcp3d_channel_layout.hpp"

using namespace std;

void MChannelLayoutTest::SetUp(const string &channel_root_name, bool channel_root_name_is_full_path)
{
    channel_root_dir_ = channel_root_name_is_full_path ? channel_root_name : mcp3d::JoinPath(test_channel_layout_dir(), channel_root_name);
    channel_layout_ = make_unique<mcp3d::MChannelLayout>(channel_root_dir_);
}

void MChannelLayoutTest::ConfigurePyrLevelDirs(const std::vector<int> &pyr_levels)
{
    ASSERT_TRUE(mcp3d::RemovePath(channel_root_dir_));
    ASSERT_TRUE(mcp3d::MakeDirectories(channel_root_dir_));
    sort(pyr_levels.begin(), pyr_levels.end());
    // level 0 is in channel_root_dir_ itself. its impossible not to have it
    ASSERT_TRUE(pyr_levels[0] == 0);
    pyr_slices_fixtures_.clear();
    for (const auto& pyr_level: pyr_levels)
        RandomConfigurePyrLevelDir(pyr_level);
    ASSERT_EQ(pyr_levels.size(), pyr_slices_fixtures_.size());
    for (int pyr_level: pyr_levels)
    {
        ASSERT_NE(pyr_slices_fixtures_.end(), pyr_slices_fixtures_.find(pyr_level));
        ASSERT_EQ(mcp3d::JoinPath(channel_root_dir_, mcp3d::MChannelLayout::PyrLevelDirName(pyr_level)), pyr_slices_fixtures_.at(pyr_level).channel_pyr_dir_);
    }
}

void MChannelLayoutTest::RandomConfigurePyrLevelDir(int pyr_level)
{
    // 50% chance for each axis to be flat.
    // iterate through slices. for each slice dir with random chance add to the slide
    // dir a file from {1.ome.tif, 1.hdf5, 1.fake}. once added stop iteration
    double p_axis_flat = 0.5, p_add_file = 0.25;
    default_random_engine generator;
    bernoulli_distribution axis_flat(p_axis_flat), add_file(p_add_file);
    uniform_int_distribution<int> number(2, 100), dim(7, 12), file_id(0, 2);
    vector<string> file_names({"1.ome.tif", "1.hdf5", "1.tif", "1.fake"});

    string pyr_dir_name = mcp3d::MChannelLayout::PyrLevelDirName(pyr_level);
    vector<int> slice_numbers, slice_dims;
    for (int i = 0; i < 3; ++i)
    {
        if (axis_flat(generator))
        {
            slice_numbers.push_back(1);
            slice_dims.push_back(mcp3d::MAX_AXIS_SLICE_DIM);
        }
        else
        {
            slice_numbers.push_back(number(generator));
            // dimension are powers of two between [128, 4096]
            slice_dims.push_back((int)mcp3d::IntPow(2, dim(generator)));
        }
    }
    pyr_slices_fixtures_.emplace(pyr_level, MChannelPyrSlicesTest{});
    pyr_slices_fixtures_.at(pyr_level).Configure(mcp3d::JoinPath(channel_root_dir_, pyr_dir_name), slice_numbers, slice_dims, true);
    bool add = false;
    for (int zid = 0; zid < slice_numbers[0], !add; ++zid)
        for (int yid = 0; yid < slice_numbers[1], !add; ++yid)
            for (int xid = 0; xid < slice_numbers[2], !add; ++xid)
            {
                add = add_file(generator);
                if (add)
                    pyr_slices_fixtures_.at(pyr_level).AddFiles(vector<string>({file_names[file_id(generator)]}), zid, yid, xid);
            }
}

void MChannelLayoutTest::ConfigureImarisDir()
{
    // flat channel root dir
    pyr_slices_fixtures_.emplace(0, MChannelPyrSlicesTest{});
    pyr_slices_fixtures_.at(0).Configure(channel_root_dir_, vector<int>(3, 1), vector<int>(3, mcp3d::MAX_AXIS_SLICE_DIM), true);
    // place symlink to imaris file
    string imaris_target = mcp3d::JoinPath(mcp3d::test_data_dir(),
                                         "src/image_formats/hdf5_format/full_sagittal_section_300uM_cropped_challenging_FusionStitcher.ims"),
           link = mcp3d::JoinPath(test_channel_layout_dir(), "1.FusionStitcher.ims");
    ASSERT_TRUE(mcp3d::Hdf5Handle(imaris_target));
    string command = mcp3d::JoinVector(vector<string>({"ln", "-s", imaris_target, link}), " ");
    mcp3d::SysCmdResult(command.c_str());
    ASSERT_TRUE(mcp3d::Hdf5Handle(link));
    // update pyr_slices_fixtures_
    pyr_slices_fixtures_.at(0).pyr_slices_.Refresh();
    int n = mcp3d::ImarisResolutionNumberOfChannels(imaris_target, 0);
    for (int i = 1; i < n; ++i)
        pyr_slices_fixtures_.emplace(n, MChannelPyrSlicesTest{pyr_slices_fixtures_.at(0)});
}

TEST(MChannelLayout, IsPyrLevelDir)
{
    EXPECT_FALSE(mcp3d::MChannelLayout::IsPyrLevelDir(""));
    EXPECT_FALSE(mcp3d::MChannelLayout::IsPyrLevelDir("pyr_level_00"));
    EXPECT_TRUE(mcp3d::MChannelLayout::IsPyrLevelDir("pyr_level_02"));
    EXPECT_TRUE(mcp3d::MChannelLayout::IsPyrLevelDir("/pyr_level_51"));
    EXPECT_TRUE(mcp3d::MChannelLayout::IsPyrLevelDir("pyr_level_01/"));
    EXPECT_TRUE(mcp3d::MChannelLayout::IsPyrLevelDir("/pyr_level_01/"));
    EXPECT_TRUE(mcp3d::MChannelLayout::IsPyrLevelDir("foo/bar//pyr_level_01/"));
    EXPECT_FALSE(mcp3d::MChannelLayout::IsPyrLevelDir("pyr_level_"));
    EXPECT_FALSE(mcp3d::MChannelLayout::IsPyrLevelDir("pyr_level_1"));
    EXPECT_FALSE(mcp3d::MChannelLayout::IsPyrLevelDir("pyr_level_111"));
    EXPECT_FALSE(mcp3d::MChannelLayout::IsPyrLevelDir("apyr_level_11"));
    EXPECT_FALSE(mcp3d::MChannelLayout::IsPyrLevelDir("pyr_level_11_"));
}

TEST(MChannelLayout, DirPyrLevel)
{
    EXPECT_EQ(0, mcp3d::MChannelLayout::DirPyrLevel("foo/pyr_level_0"));
    EXPECT_EQ(7, mcp3d::MChannelLayout::DirPyrLevel("foo/pyr_level_07"));
    EXPECT_EQ(7, mcp3d::MChannelLayout::DirPyrLevel("foo/pyr_level_07/"));
    EXPECT_EQ(77, mcp3d::MChannelLayout::DirPyrLevel("foo/pyr_level_77"));
    EXPECT_EQ(0, mcp3d::MChannelLayout::DirPyrLevel("foo/pyr_level_777"));
    EXPECT_EQ(0, mcp3d::MChannelLayout::DirPyrLevel("foo/pyr_level_77/a"));
}

TEST(MChannelLayout, PyrLevelDirName)
{
    EXPECT_THROW(mcp3d::MChannelLayout::PyrLevelDirName(100), mcp3d::MCPAssertionError);
    EXPECT_EQ("", mcp3d::MChannelLayout::PyrLevelDirName(0));
    EXPECT_EQ("pyr_level_02", mcp3d::MChannelLayout::PyrLevelDirName(2));
    EXPECT_EQ("pyr_level_90", mcp3d::MChannelLayout::PyrLevelDirName(90));
}

TEST_F(MChannelLayoutTest, ReadChannelLayout)
{
    auto DrawLevels = [](){
        unordered_set<int> levels({0});
        default_random_engine generator;
        uniform_int_distribution<int> level_range(1, 99);
        while (levels.size() < 10)
            levels.insert(level_range(generator));
        return levels;
    };

    unordered_set<int> seen_levels, current_levels;
    // randomly generate non imaris pyr level directories and configure them
    for (int i = 0; i < 10; ++i)
    {
        current_levels = DrawLevels();
        seen_levels.insert(current_levels.begin(), current_levels.end());
        ConfigurePyrLevelDirs(vector<int>(current_levels.begin(), current_levels.end()));
        for (int pyr_level: seen_levels)
        {
            if (current_levels.find(pyr_level) == current_levels.end())
            {
                ASSERT_FALSE(channel_layout_->HasPyrLevelDir(pyr_level));
                ASSERT_FALSE(channel_layout_->HasPyrLevel(pyr_level));
                ASSERT_EQ(pyr_slices_fixtures_.at(pyr_level).pyr_slices(), channel_layout_->channel_pyr_slices(pyr_level));
            }
        }

    }
}

/*
TEST(MVolumeLayout, Constructor)
{
    string imaris_data_dir = mcp3d::JoinPath(mcp3d::test_data_dir(), "image_formats/hdf5_format");
    string tiff_data_dir = mcp3d::JoinPath(mcp3d::test_data_dir(), "image_formats/tiff_format/large_files_read");
    vector<string> channel_dir_names({"ch0", "ch1"});

    // scope type unknown, file format found by mcp3d::SearchFileFormat
    mcp3d::MVolumeLayout tiff_data_dir_layout(tiff_data_dir);
    EXPECT_EQ(tiff_data_dir, tiff_data_dir_layout.volume_root_dir());
    EXPECT_EQ(mcp3d::FileFormat::TIFF, tiff_data_dir_layout.file_format());
    EXPECT_EQ(vector<string>({""}), tiff_data_dir_layout.channel_dir_names());
    EXPECT_EQ(1, tiff_data_dir_layout.n_channels());
    EXPECT_EQ(mcp3d::ScopeType::UNKNOWN, tiff_data_dir_layout.scope_type());
    EXPECT_EQ(mcp3d::TissueType::UNKNOWN, tiff_data_dir_layout.tissue_type());
    EXPECT_EQ(mcp3d::CuttingPlane::UNKNOWN, tiff_data_dir_layout.cutting_plane());
    EXPECT_EQ(-1, tiff_data_dir_layout.magnification());
    EXPECT_EQ(-1, tiff_data_dir_layout.thickness());
    EXPECT_EQ("", tiff_data_dir_layout.slice_name());
    EXPECT_EQ("", tiff_data_dir_layout.slide_name());
    EXPECT_EQ("", tiff_data_dir_layout.region_name());

    // scope type unknown, file format found by mcp3d::SearchFileFormat
    mcp3d::MVolumeLayout imaris_data_dir_layout(imaris_data_dir);
    EXPECT_EQ(imaris_data_dir, imaris_data_dir_layout.volume_root_dir());
    EXPECT_EQ(mcp3d::FileFormat::IMARIS, imaris_data_dir_layout.file_format());
    EXPECT_EQ(vector<string>({""}), imaris_data_dir_layout.channel_dir_names());
    EXPECT_EQ(1, imaris_data_dir_layout.n_channels());
    EXPECT_EQ(mcp3d::ScopeType::UNKNOWN, imaris_data_dir_layout.scope_type());
    EXPECT_EQ(mcp3d::TissueType::UNKNOWN, imaris_data_dir_layout.tissue_type());
    EXPECT_EQ(mcp3d::CuttingPlane::UNKNOWN, imaris_data_dir_layout.cutting_plane());
    EXPECT_EQ(-1, imaris_data_dir_layout.magnification());
    EXPECT_EQ(-1, imaris_data_dir_layout.thickness());
    EXPECT_EQ("", imaris_data_dir_layout.slice_name());
    EXPECT_EQ("", imaris_data_dir_layout.slide_name());
    EXPECT_EQ("", imaris_data_dir_layout.region_name());

    mcp3d::MakeDirectories(test_volume_layout_dir());

    // dragonfly directory, all fields should be parsed
    // dragonfly_dir0 -> imaris_data_dir
    string dragonfly_dir0 = mcp3d::JoinPath(test_volume_layout_dir(), "Raw/Dragonfly_10x/WholeBrain/1_05_400um_BLA_coronal/a");
    mcp3d::MakeDirectories(mcp3d::ParentDir(dragonfly_dir0));
    mcp3d::SysCmdResult("ln -s -L " + imaris_data_dir + " " + dragonfly_dir0);
    mcp3d::MVolumeLayout dragfly_dir0_layout(dragonfly_dir0);
    EXPECT_EQ(dragonfly_dir0, dragfly_dir0_layout.volume_root_dir());
    EXPECT_EQ(mcp3d::FileFormat::IMARIS, dragfly_dir0_layout.file_format());
    EXPECT_EQ(vector<string>({""}), dragfly_dir0_layout.channel_dir_names());
    EXPECT_EQ(1, dragfly_dir0_layout.n_channels());
    EXPECT_EQ(mcp3d::ScopeType::DRAGONFLY, dragfly_dir0_layout.scope_type());
    EXPECT_EQ(mcp3d::TissueType::WHOLE_BRAIN, dragfly_dir0_layout.tissue_type());
    EXPECT_EQ(mcp3d::CuttingPlane::CORONAL, dragfly_dir0_layout.cutting_plane());
    EXPECT_EQ(10, dragfly_dir0_layout.magnification());
    EXPECT_EQ(400, dragfly_dir0_layout.thickness());
    EXPECT_EQ("1_05_400um_BLA_coronal", dragfly_dir0_layout.slice_name());
    EXPECT_EQ("1_05", dragfly_dir0_layout.slide_name());
    EXPECT_EQ("BLA", dragfly_dir0_layout.region_name());

    // dragonfly file path, all fields should be parsed
    string dragonfly_path0 = mcp3d::JoinPath(dragonfly_dir0, "full_sagittal_section_300uM_cropped_challenging_FusionStitcher.ims");
    mcp3d::MVolumeLayout dragfly_path0_layout(dragonfly_path0);
    EXPECT_EQ(dragonfly_dir0, dragfly_path0_layout.volume_root_dir());
    EXPECT_EQ(mcp3d::FileFormat::IMARIS, dragfly_path0_layout.file_format());
    EXPECT_EQ(vector<string>({""}), dragfly_path0_layout.channel_dir_names());
    EXPECT_EQ(1, dragfly_path0_layout.n_channels());
    EXPECT_EQ(mcp3d::ScopeType::DRAGONFLY, dragfly_path0_layout.scope_type());
    EXPECT_EQ(mcp3d::TissueType::WHOLE_BRAIN, dragfly_path0_layout.tissue_type());
    EXPECT_EQ(mcp3d::CuttingPlane::CORONAL, dragfly_path0_layout.cutting_plane());
    EXPECT_EQ(10, dragfly_path0_layout.magnification());
    EXPECT_EQ(400, dragfly_path0_layout.thickness());
    EXPECT_EQ("1_05_400um_BLA_coronal", dragfly_path0_layout.slice_name());
    EXPECT_EQ("1_05", dragfly_path0_layout.slide_name());
    EXPECT_EQ("BLA", dragfly_path0_layout.region_name());

    // dragonfly directory, only some fields can be parsed. give channel_dir_names,
    // which should be ignored when imaris file format found
    // dragonfly_dir1 -> imaris_data_dir
    string dragonfly_dir1 = mcp3d::JoinPath(test_volume_layout_dir(), "Raw/Dragonfly_10x/1_05");
    mcp3d::MakeDirectories(mcp3d::ParentDir(dragonfly_dir1));
    mcp3d::SysCmdResult("ln -s -L " + imaris_data_dir + " " + dragonfly_dir1);
    mcp3d::MVolumeLayout dragfly_dir1_layout(dragonfly_dir1, channel_dir_names);
    EXPECT_EQ(dragonfly_dir1, dragfly_dir1_layout.volume_root_dir());
    EXPECT_EQ(mcp3d::FileFormat::IMARIS, dragfly_dir1_layout.file_format());
    EXPECT_EQ(vector<string>({""}), dragfly_dir1_layout.channel_dir_names());
    EXPECT_EQ(1, dragfly_dir1_layout.n_channels());
    EXPECT_EQ(mcp3d::ScopeType::DRAGONFLY, dragfly_dir1_layout.scope_type());
    EXPECT_EQ(mcp3d::TissueType::UNKNOWN, dragfly_dir1_layout.tissue_type());
    EXPECT_EQ(mcp3d::CuttingPlane::UNKNOWN, dragfly_dir1_layout.cutting_plane());
    EXPECT_EQ(10, dragfly_dir1_layout.magnification());
    EXPECT_EQ(-1, dragfly_dir1_layout.thickness());
    EXPECT_EQ("", dragfly_dir1_layout.slice_name());
    EXPECT_EQ("", dragfly_dir1_layout.slide_name());
    EXPECT_EQ("", dragfly_dir1_layout.region_name());

    // lightsheet data directory ending in lightsheet channel dir name. no files
    // under the directory.
    // give channel_dir_names, which should be ignored
    string lightsheet_dir0 = mcp3d::JoinPath(test_volume_layout_dir(), "Raw/Lightsheet_10x/SpinalCord/Ex_0_Em_0");
    mcp3d::MakeDirectories(lightsheet_dir0);
    mcp3d::MakeDirectories(mcp3d::JoinPath(test_volume_layout_dir(), "Raw/Lightsheet_10x/SpinalCord/Ex_1_Em_1"));
    mcp3d::MVolumeLayout lightsheet_dir0_layout(mcp3d::ParentDir(lightsheet_dir0), channel_dir_names);
    EXPECT_EQ(mcp3d::ParentDir(lightsheet_dir0),
              lightsheet_dir0_layout.volume_root_dir());
    EXPECT_EQ(mcp3d::FileFormat::UNKNOWN, lightsheet_dir0_layout.file_format());
    EXPECT_EQ(vector<string>({"Ex_0_Em_0", "Ex_1_Em_1"}), lightsheet_dir0_layout.channel_dir_names());
    EXPECT_EQ(2, lightsheet_dir0_layout.n_channels());
    EXPECT_EQ(mcp3d::ScopeType::LIGHTSHEET, lightsheet_dir0_layout.scope_type());
    EXPECT_EQ(mcp3d::TissueType::SPINAL_CORD, lightsheet_dir0_layout.tissue_type());
    EXPECT_EQ(mcp3d::CuttingPlane::UNKNOWN, lightsheet_dir0_layout.cutting_plane());
    EXPECT_EQ(10, lightsheet_dir0_layout.magnification());
    EXPECT_EQ(-1, lightsheet_dir0_layout.thickness());
    EXPECT_EQ("", lightsheet_dir0_layout.slice_name());
    EXPECT_EQ("", lightsheet_dir0_layout.slide_name());
    EXPECT_EQ("", lightsheet_dir0_layout.region_name());

    // lightsheet with no scope magnification combination, resulting in unknown scope type
    // lightsheet_dir1 -> tiff_data_dir
    string lightsheet_dir1 = mcp3d::JoinPath(test_volume_layout_dir(), "Raw/Lightsheet/SpinalCord/Ex_0_Em_0");
    mcp3d::MakeDirectories(mcp3d::ParentDir(lightsheet_dir1));
    mcp3d::SysCmdResult("ln -s -L " + tiff_data_dir + " " + lightsheet_dir1);
    mcp3d::MVolumeLayout lightsheet_dir1_layout(lightsheet_dir1);
    EXPECT_EQ(lightsheet_dir1, lightsheet_dir1_layout.volume_root_dir());
    EXPECT_EQ(mcp3d::FileFormat::TIFF, lightsheet_dir1_layout.file_format());
    EXPECT_EQ(vector<string>({""}), lightsheet_dir1_layout.channel_dir_names());
    EXPECT_EQ(1, lightsheet_dir1_layout.n_channels());
    EXPECT_EQ(mcp3d::ScopeType::UNKNOWN, lightsheet_dir1_layout.scope_type());
    EXPECT_EQ(mcp3d::TissueType::UNKNOWN, lightsheet_dir1_layout.tissue_type());
    EXPECT_EQ(mcp3d::CuttingPlane::UNKNOWN, lightsheet_dir1_layout.cutting_plane());
    EXPECT_EQ(-1, lightsheet_dir1_layout.magnification());
    EXPECT_EQ(-1, lightsheet_dir1_layout.thickness());
    EXPECT_EQ("", lightsheet_dir1_layout.slice_name());
    EXPECT_EQ("", lightsheet_dir1_layout.slide_name());
    EXPECT_EQ("", lightsheet_dir1_layout.region_name());

    // remove symlink to tiff_data_dir, and make a normal directory
    mcp3d::RemovePath(lightsheet_dir1, false);  // do not remove link target
    mcp3d::MakeDirectories(lightsheet_dir1);
    lightsheet_dir1_layout = mcp3d::MVolumeLayout(lightsheet_dir1);
    EXPECT_EQ(lightsheet_dir1, lightsheet_dir1_layout.volume_root_dir());
    EXPECT_EQ(mcp3d::FileFormat::UNKNOWN, lightsheet_dir1_layout.file_format());
    EXPECT_EQ(vector<string>({""}), lightsheet_dir1_layout.channel_dir_names());
    EXPECT_EQ(1, lightsheet_dir1_layout.n_channels());
    EXPECT_EQ(mcp3d::ScopeType::UNKNOWN, lightsheet_dir1_layout.scope_type());
    EXPECT_EQ(mcp3d::TissueType::UNKNOWN, lightsheet_dir1_layout.tissue_type());
    EXPECT_EQ(mcp3d::CuttingPlane::UNKNOWN, lightsheet_dir1_layout.cutting_plane());
    EXPECT_EQ(-1, lightsheet_dir1_layout.magnification());
    EXPECT_EQ(-1, lightsheet_dir1_layout.thickness());
    EXPECT_EQ("", lightsheet_dir1_layout.slice_name());
    EXPECT_EQ("", lightsheet_dir1_layout.slide_name());
    EXPECT_EQ("", lightsheet_dir1_layout.region_name());

    // unknown scope type, channel directories do not exist
    lightsheet_dir1_layout = mcp3d::MVolumeLayout(lightsheet_dir1, channel_dir_names);
    EXPECT_EQ(lightsheet_dir1, lightsheet_dir1_layout.volume_root_dir());
    EXPECT_EQ(mcp3d::FileFormat::UNKNOWN, lightsheet_dir1_layout.file_format());
    EXPECT_EQ(vector<string>({""}), lightsheet_dir1_layout.channel_dir_names());
    EXPECT_EQ(1, lightsheet_dir1_layout.n_channels());
    EXPECT_EQ(mcp3d::ScopeType::UNKNOWN, lightsheet_dir1_layout.scope_type());
    EXPECT_EQ(mcp3d::TissueType::UNKNOWN, lightsheet_dir1_layout.tissue_type());
    EXPECT_EQ(mcp3d::CuttingPlane::UNKNOWN, lightsheet_dir1_layout.cutting_plane());
    EXPECT_EQ(-1, lightsheet_dir1_layout.magnification());
    EXPECT_EQ(-1, lightsheet_dir1_layout.thickness());
    EXPECT_EQ("", lightsheet_dir1_layout.slice_name());
    EXPECT_EQ("", lightsheet_dir1_layout.slide_name());
    EXPECT_EQ("", lightsheet_dir1_layout.region_name());

    // unknown scope type, channel directories exist, but no files under them
    for (const auto& channel_dir_name: channel_dir_names)
        mcp3d::MakeDirectories(mcp3d::JoinPath({lightsheet_dir1, channel_dir_name}));
    lightsheet_dir1_layout = mcp3d::MVolumeLayout(lightsheet_dir1, channel_dir_names);
    EXPECT_EQ(lightsheet_dir1, lightsheet_dir1_layout.volume_root_dir());
    EXPECT_EQ(mcp3d::FileFormat::UNKNOWN, lightsheet_dir1_layout.file_format());
    EXPECT_EQ(channel_dir_names, lightsheet_dir1_layout.channel_dir_names());
    EXPECT_EQ(2, lightsheet_dir1_layout.n_channels());
    EXPECT_EQ(mcp3d::ScopeType::UNKNOWN, lightsheet_dir1_layout.scope_type());
    EXPECT_EQ(mcp3d::TissueType::UNKNOWN, lightsheet_dir1_layout.tissue_type());
    EXPECT_EQ(mcp3d::CuttingPlane::UNKNOWN, lightsheet_dir1_layout.cutting_plane());
    EXPECT_EQ(-1, lightsheet_dir1_layout.magnification());
    EXPECT_EQ(-1, lightsheet_dir1_layout.thickness());
    EXPECT_EQ("", lightsheet_dir1_layout.slice_name());
    EXPECT_EQ("", lightsheet_dir1_layout.slide_name());
    EXPECT_EQ("", lightsheet_dir1_layout.region_name());

    mcp3d::RemovePath(test_volume_layout_dir());
}
 */