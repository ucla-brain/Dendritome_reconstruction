//
// Created by muyezhu on 9/28/20.
//

#ifndef MCP3D_TEST_MCP3D_CHANNEL_LAYOUT_HPP
#define MCP3D_TEST_MCP3D_CHANNEL_LAYOUT_HPP

#include <memory>
#include <vector>
#include <unordered_map>
#include "common/mcp3d_common.hpp"
#include "image_layout/mcp3d_channel_layout.hpp"
#include "test_mcp3d_channel_pyr_slices.hpp"

// test fixture. will use ASSERT_* to abort invalid testing conditions
class MChannelLayoutTest: public ::testing::Test
{
protected:
    MChannelLayoutTest(): channel_root_dir_(std::string{}), channel_layout_(nullptr), pyr_slices_fixtures_(std::unordered_map<int, MChannelPyrSlicesTest>{}) {};

    void SetUp(const std::string& channel_root_name, bool channel_root_name_is_full_path = false);

    /// remove channel_root_dir_ and repopulate it. clear pyr_slices_fixtures_ and repopulate from newly generated channel_root_dir_ layout.
    /// sort pyr_levels: MChannelPyrSlicesTest class will remove its channel_pyr_dir in Configure function. if pyr level 0 is not configured first,
    /// it will cause earlier configured level directories to be removed.
    /// will not create stitched imaris directories
    void ConfigurePyrLevelDirs(const std::vector<int> &pyr_levels);

    /// create stitched imaris directory
    void ConfigureImarisDir();

    std::string test_channel_layout_dir()
    { return mcp3d::JoinPath(mcp3d::test_data_dir(), "image_layout", "channel_layout"); }

    std::string channel_root_dir_;
    /// using a pointer because fixtures can not be passed constructor arguments in TEST_F
    /// MChannelLayout class has deleted default constructor and can not be instantiated by the default MChannelLayoutTest constructor
    std::unique_ptr<mcp3d::MChannelLayout> channel_layout_;
    std::unordered_map<int, MChannelPyrSlicesTest> pyr_slices_fixtures_;

private:
    void RandomConfigurePyrLevelDir(int pyr_level);

};

#endif //MCP3D_TEST_MCP3D_CHANNEL_LAYOUT_HPP
