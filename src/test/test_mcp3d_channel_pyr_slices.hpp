//
// Created by muyezhu on 9/28/20.
//

#ifndef MCP3D_TEST_MCP3D_CHANNEL_PYR_SLICES_HPP
#define MCP3D_TEST_MCP3D_CHANNEL_PYR_SLICES_HPP

#include <string>
#include <vector>
#include <gtest/gtest.h>

// test fixture. will use ASSERT_* to abort invalid testing conditions
class MChannelPyrSlicesTest: public ::testing::Test
{
public:
    /// if !dir_name_is_full_path, set channel_pyr_dir to test_channel_pyr_slices_dir()/dir_name. otherwise set channel_pyr_dir to dir_name
    /// Configure, AddDirs, AddFiles, Removedir and RemoveFiles updates pyr_slices_ to reflect file/directory changes
    void Configure(const std::string &dir_name, const std::vector<int> &slice_numbers, const std::vector<int> &slice_dims, bool dir_name_is_full_path = false);

    void AddDirs(const std::vector<std::string> &dir_paths);

    void AddFiles(const std::vector<std::string>& file_names, int slice_zid, int slice_yid, int slice_xid);

    // if any slice id is negative, delete all dirs and files identified by the previous non negative slice ids
    void RemoveDir(int slice_zid, int slice_yid, int slice_xid);

    void RemoveFiles(const std::vector<std::string>& file_names, int slice_zid, int slice_yid, int slice_xid);

    void RefreshPyrSlices();

    const mcp3d::MChannelPyrSlices& pyr_slices() const
    { return pyr_slices_; }

protected:
    /// this function exists because slice dir name is required before pyr_slices_ starts to read the slices
    void GetSliceDir(int slice_zid, int slice_yid, int slice_xid, std::string& slice_dir) const;

    void TearDown() override
    { CleanUp(); }

    std::string test_channel_pyr_slices_dir()
    { return mcp3d::JoinPath(mcp3d::test_data_dir(), "image_layout", "channel_pyr_slices"); }

    std::string channel_pyr_dir_;
    std::vector<int> slice_dims_, slice_numbers_;
    mcp3d::MChannelPyrSlices pyr_slices_;

private:
    void CleanUp()
    { mcp3d::RemovePath(channel_pyr_dir_); }
};

#endif //MCP3D_TEST_MCP3D_CHANNEL_PYR_SLICES_HPP
