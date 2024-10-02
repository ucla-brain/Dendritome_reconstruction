//
// Created by muyezhu on 6/17/19.
//

#ifndef MCP3D_MCP3D_IMAGE_INFO_HPP
#define MCP3D_MCP3D_IMAGE_INFO_HPP

#include <utility>
#include <unordered_set>
#include <unordered_map>
#include "image_layout/mcp3d_volume_layout.hpp"
#include "image_interface/mcp3d_voxel_types.hpp"
#include "mcp3d_channel_info.hpp"

namespace mcp3d
{

class MImageInfo
{
public:
    MImageInfo() = default;

    /// AssertValid() before constructor exists
    explicit MImageInfo(const std::string& volume_path, const std::vector<std::string>& channel_dir_names = std::vector<std::string>{});

    /// transfer ownership of channel_layouts under volume_layout to volume_layout_ member through MVolumeLayout move constructor
    explicit MImageInfo(MVolumeLayout&& volume_layout): volume_layout_(volume_layout), channel_infos_(std::unordered_map<std::string, MChannelInfo>{})
    { AssertValid(); };

    /// deep copy of volume_layout_. no ownership sharing between this and other
    MImageInfo(const MImageInfo& other): volume_layout_(other.volume_layout_), channel_infos_(other.channel_infos_)
    { AssertValid(); };

    bool operator==(const MImageInfo& other) const
    { return volume_layout_ == other.volume_layout_ && channel_infos_ == other.channel_infos_; }

    /// return true if volume_layout_.empty()
    bool empty() const
    { return volume_layout_.empty(); }

    /// channel_pyr_info will be moved from.
    /// AssertValid()
    void AddChannelPyrInfo(const std::string& channel_name, MChannelPyrInfo&& channel_pyr_info);

    /// channel_pyr_info will be moved from
    void AddChannelPyrInfo(int channel_number, MChannelPyrInfo&& channel_pyr_info)
    { AddChannelPyrInfo(volume_layout_.channel_name(channel_number, 0), std::move(channel_pyr_info)); }

    /// AssertValid()
    /// call MChannelInfo::Save for all entries in channel_infos_
    /// all channel info json files are named __channel_info__.json
    /// for imaris formats where channel dir names are empty string, channel info
    /// json for all channels will have the same path, therefore only channel 0
    /// needs to write its __channel_info__.json to disk
    void Save();

    void ClearChannelPyrInfos(const std::string& channel_name)
    { channel_info(channel_name).ClearChannelPyrInfos(); }

    void ClearChannelPyrInfos(int channel_number, int time = 0)
    { channel_info(channel_number, time).ClearChannelPyrInfos(); }

    /// if volume_layout_.HasChannelName(channel_name) is false, return false
    /// else return channel_infos_.at(channel_name).HasPyrInfo(pyr_level)
    bool HasChannelPyrInfo(const std::string& channel_name, int pyr_level) const;

    bool HasChannelPyrInfo(int channel_number, int pyr_level) const;

    int n_pyr_levels(const std::string& channel_name) const;

    int n_pyr_levels(int channel_number = 0, int time = 0) const
    { return n_pyr_levels(volume_layout_.channel_name(channel_number, time)); }

    int n_pyr_infos(const std::string& channel_name) const;

    int n_pyr_infos(int channel_number = 0, int time = 0) const
    { return n_pyr_infos(volume_layout_.channel_name(channel_number, time)); }

    /// returns the highest (coarsest) pyramid level that can be mapped correctly to level 0 image (based on image xy dimensions)
    /// implementation assumes that any MChannelPyrInfo represents pyr level volumes mappable to level 0
    /// z dimension is allowed to maintain its pyr ratio at any given level, while xy dimensions must be halved from level to level + 1
    /// if no MChannelPyrInfo has been read across any channels, return -1 since the actual result can not be known
    int MaxMappablePyrLevel() const;

    int zdim(int pyr_level = 0) const
    { return xyz_dims(pyr_level)[0]; }

    int ydim(int pyr_level = 0) const
    { return xyz_dims(pyr_level)[1]; }

    int xdim(int pyr_level = 0) const
    { return xyz_dims(pyr_level)[2]; }

    /// return xyz_dims(pyr_level) from any MChannelInfo that has pyr_level in channel_pyr_infos_
    /// if no MChannelInfo has pyr_level in channel_pyr_infos_, return {0, 0, 0}
    std::vector<int> xyz_dims(int pyr_level = 0) const;

    /// return dimensions of time, channel, z, y, x = {1, volume_layout_.n_channels(), xyz_dims(0)}
    /// if not level MChannelPyrInfo is present in any MChannelInfo, the zyx dimensions are returned as 0s
    std::vector<int> dims() const
    { return std::vector<int>({1, volume_layout_.n_channels(), zdim(0), ydim(0), xdim(0)}); }

    // getters
    const MVolumeLayout& volume_layout() const
    { return volume_layout_; }

    MVolumeLayout& volume_layout()
    { return volume_layout_; }

    const MChannelInfo& channel_info(const std::string& channel_name) const;

    MChannelInfo& channel_info(const std::string& channel_name);

    const MChannelInfo& channel_info(int channel_number = 0, int time = 0) const
    { return channel_info(volume_layout_.channel_name(channel_number, time)); }

    MChannelInfo& channel_info(int channel_number = 0, int time = 0)
    { return channel_info(volume_layout_.channel_name(channel_number, time)); }

    const MChannelPyrInfo& channel_pyr_info(const std::string& channel_name, int pyr_level) const
    { return channel_info(channel_name).channel_pyr_info(pyr_level); }

    const MChannelPyrInfo& channel_pyr_info(int channel_number, int pyr_level) const
    { return channel_info(channel_number, 0).channel_pyr_info(pyr_level); }

private:
    /// return if empty()
    /// (1) volume_layout_.AssertValid
    /// (2) channel_infos_ has volume_layout_.channel_names_ as keys and nothing else
    ///     channel_infos[channel_name].channel_layout_ store the same pointer as volume_layout_.channel_layouts_[channel_name]
    ///     call channel_infos[channel_name].AssertValid
    /// (3) all MChannelPyrInfo has same xyz_dims at same pyr levels (other MChannelPyrInfo members such as voxel type, MChannelPyrSlices,
    ///     chunk xyz dimensions, file format (if not imaris) can be different)
    void AssertValid() const;

    /// construct MChannelInfo with volume_layout_.channel_layouts_'s shared_ptr<MChannelLayout> instances
    /// the constructed MChannelInfo share ownership of MChannelLayout instance with volume_layout_ and
    /// has empty channel_pyr_infos_
    void InitChannelInfos();

    MVolumeLayout volume_layout_;
    /// channel name: MChannelInfo.
    /// MChannelInfo instances share ownership of MChannelLayout instances with volume_layout_
    /// channel_infos_ should not be cleared. since volume_layout_ has static channel_names_
    /// once constructed, conceptually there should always be the same set of channel_name: MChannelInfo
    /// pairs associated with the image volume. if the files/directories under the image volume
    /// changes, these changes should be reflected through MVolumeLayout::RefreshChannelLayouts
    /// as well as image format specific MImageIO operations to read MChannelPyrInfo instances
    std::unordered_map<std::string, MChannelInfo> channel_infos_;
};

}

#endif //MCP3D_MCP3D_IMAGE_INFO_HPP
