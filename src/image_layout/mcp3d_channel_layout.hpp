//
// Created by muyezhu on 3/23/18.
//

#ifndef MCP3D_MCP3D_IMAGE_LAYOUT_HPP
#define MCP3D_MCP3D_IMAGE_LAYOUT_HPP

#include <regex>
#include <unordered_map>
#include "common/mcp3d_common.hpp"
#include "image_interface/mcp3d_file_formats.hpp"
#include "mcp3d_channel_pyr_slices.hpp"

namespace mcp3d
{

/// member channel_root_dir_ is const and can not be modified once constructed
class MChannelLayout
{
public:
    MChannelLayout() = delete;

    /// call ReadChannelLayout
    explicit MChannelLayout(const std::string& channel_root_dir);

    /// deep copy of MChannelPyrSlices instances managed by other.channel_pyr_slices_
    /// other's shared_ptr do not share object management with this instance
    MChannelLayout(const MChannelLayout& other);

    /// move construction. pyr_channel_slices_ will share object management with
    /// shared_ptr that previously shares object management with other.
    /// other.pyr_channel_slices_ will be empty
    MChannelLayout(MChannelLayout&& other) noexcept: channel_root_dir_(other.channel_root_dir_), channel_pyr_slices_(other.channel_pyr_slices_) {}

    /// update the current folder structure on disk into channel_pyr_slices_
    /// call ClearPyrSlices.
    /// for imaris format, only pyr level 0 directory exists. modify channel_pyr_slices_ for
    /// other pyr levels explicitly. all MChannelPyrSlices instances have channel_pyr_dir_ = channel_root_dir_
    /// for non imaris format, at each pyr level, if the pyr level is already a key in pyr_slices_,
    /// call MChannelPyrSlices::set_channel_pyr_dir. otherwise, construct
    /// shared_ptr managing MChannelPyrSlices constructed from channel_pyr_dir.
    /// if a pyr level was in pyr_slices_ but no longer has a corresponding pyr
    /// level dir on disk (potentially removed), call Clear() on the MChannelPyrSlices instance
    void ReadChannelLayout();

    /// deep comparison of MChannelPyrSlices instances managed by shared_ptr
    /// for each pyr_level in channel_pyr_slices_, if pyr_level not in other.channel_pyr_slices_,
    /// channel_pyr_slices_[pyr_level]->empty() should be true. otherwise
    /// the two MChannelPyrSlices object should be equal. same process applies
    /// for each pyr_level in other.channel_pyr_slices_
    bool operator==(const MChannelLayout& other) const;

    bool operator!=(const MChannelLayout& other) const
    { return !(*this == other); }

    /// return true if pyr_level found in channel_pyr_slices_, and channel_pyr_slices_[pyr_level] is not empty
    /// note that no key value pairs will be erased from channel_pyr_slices_. if a pyr level dir is removed
    /// from disk this is reflected as channel_pyr_slices_[pyr_level]->Clear()
    bool HasPyrLevel(int pyr_level) const;

    /// return true if the directory channel_root_dir_/PyrLevelDirName(pyr_level) exists
    /// note that for imaris format, HasPyrLevel will be true for all levels, while HasPyrLevelDir will be true for level 0
    bool HasPyrLevelDir(int pyr_level) const
    { return mcp3d::IsDir(mcp3d::JoinPath({channel_root_dir_, MChannelLayout::PyrLevelDirName(pyr_level)})); }

    /// asserts channel_root_dir exists. level 0 always exists and is represented
    /// by string(). additionally dir names under channel_root_dir that return true from IsPyrLevelDir
    /// will be included. the output vector is sorted (string() being the first entry)
    std::vector<std::string> PyrLevelDirs(bool full_path = true) const;

    /// return size of PyrLevelDirs vector
    int NumberOfPyrLevelDirs() const
    { return (int)PyrLevelDirs(false).size(); }

    /// return sorted vector of pyr level keys of channel_pyr_slices_ for which the mapped instance is not empty
    /// this returns different number of items than PyrLevelDirs for imaris format
    std::vector<int> pyr_levels() const;

    int n_pyr_levels() const
    { return (int)pyr_levels().size(); }

    /// pyr_level should be within [0, 99]
    /// if pyr_level = 0, return "", otherwise return pyr_level_xx
    static std::string PyrLevelDirName(int pyr_level);

    /// does not assert dir exists
    /// returns true if match with pyr_level_dir_pattern() found
    static bool IsPyrLevelDir(const std::string &dir);

    /// if no match to pyr_level_dir_pattern(), return 0, as we'll assume this
    /// directory is the channel root. otherwise return the xx in pyr_level_xx.
    /// does not assert dir exists
    static int DirPyrLevel(const std::string &dir);

    FileFormat file_format(int pyr_level) const;

    std::string channel_root_dir() const
    { return channel_root_dir_; }

    /// if instance does not have pyr_level, call ReadChannelLayout to refresh
    /// channel layout. assert instance has pyr_level
    std::string channel_pyr_dir(int pyr_level);

    const MChannelPyrSlices& channel_pyr_slices(int pyr_level) const
    { return *(channel_pyr_slices_ptr(pyr_level)); }

    const std::shared_ptr<MChannelPyrSlices>& channel_pyr_slices_ptr(int pyr_level) const;

    static std::string pyr_level_dir_prefix()
    { return "pyr_level_"; }

    // matches pyr_level_xx, with the exception of pyr_level_00
    static std::regex pyr_level_dir_pattern()
    { return std::regex("^pyr_level_(?!00)([0-9]{2})$"); }

    friend class MChannelInfo;
    friend class MVolumeLayout;

private:
    /// no shared_ptr in channel_pyr_slices_ should be storing nullptr
    /// HasPyrLevel(0) is true
    /// if imaris at level 0, HasPyrLevel(pyr_level) true for all imaris resolution levels, but no
    /// pyr level dirs for levels 1 and greater should exist. all MChannelPyrSlices should have identical
    /// channel_pyr_dir_ and be of FileFormat::IMARIS.
    /// else no MChannelPyrSlices should be imaris format. HasPyrLevel and HasPyrLevelDir should return same results for all inputs
    void AssertValid() const;

    /// does not call channel_pyr_slices_.clear(). call Clear() on MChannelPyrSlices
    /// instances instead. this set MChannelPyrSlices to empty status rather than reducing use counts of the shared_ptr.
    /// if channel_pyr_slices_ is cleared instead, there can be other shared_ptr managing the
    /// MChannelPyrslices instances (e.g. MChannelPyrInfo), causing desynchronized
    /// view in MChannelInfo when pyr level directories have changed and a
    /// new ReadChannelLayout() is called. if n is not a key of channel_pyr_slices_, or
    /// channel_pyr_slices_[n]->empty() is true, pyr level n does not exist under channel_root_dir_
    void ClearChannelPyrSlices();

    const std::string channel_root_dir_;
    // resolution level: MChannelPyrSlices
    std::unordered_map<int, std::shared_ptr<MChannelPyrSlices>> channel_pyr_slices_;
};

}

#endif //MCP3D_MCP3D_IMAGE_LAYOUT_HPP
