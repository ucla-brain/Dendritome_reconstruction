//
// Created by muyezhu on 3/1/19.
//

#ifndef MCP3D_MCP3D_PROJECT_LAYOUT_HPP
#define MCP3D_MCP3D_PROJECT_LAYOUT_HPP

#include <memory>
#include "image_interface/mcp3d_voxel_types.hpp"
#include "image_interface/mcp3d_file_formats.hpp"
#include "mcp3d_channel_layout.hpp"

namespace mcp3d
{
/* a single imaris file includes all channels as well as all resolution levels
 * a directory of tiff images include a single channel at a single resolution level
 */
/// volume_root_, channel_names_ and channel_dir_names_ can not be changed once
/// instance construction finishes
class MVolumeLayout
{
public:
    MVolumeLayout(): volume_root_dir_(std::string{}), channel_dir_names_(std::vector<std::string>{}), channel_names_(std::vector<std::string>{}),
                     channel_layouts_(std::unordered_map<std::string, std::shared_ptr<MChannelLayout>>{}), constructed_(true) {};

    /// if volume_path is not a directory, throw exception
    /// calls SetChannelDirNames, SetChannelFileFormats
    explicit MVolumeLayout(const std::string &volume_root_dir, const std::vector<std::string>& channel_dir_names = std::vector<std::string>{});

    /// copy constructor (deep copy of channel_layouts_, do not shared ownership of MChannelLayout instances)
    MVolumeLayout(const MVolumeLayout& other);

    /// move constructor. move assign other.channel_layouts_ to channel_layouts_
    MVolumeLayout(MVolumeLayout&& other) noexcept:
            volume_root_dir_(other.volume_root_dir_), channel_dir_names_(other.channel_dir_names_), channel_names_(other.channel_names_),
            channel_layouts_(other.channel_layouts_), constructed_(other.constructed_)
    { MCP3D_ASSERT(constructed_) }


    /// return a string representation of the instance
    void Print() const;

    /// deep comparison of channel_layout_
    bool operator==(const MVolumeLayout& other) const;

    bool operator!=(const MVolumeLayout& other) const
    { return !(*this == other); }

    bool empty() const
    { return channel_names_.empty(); }

    /// return true if channel_name is found in channel_names_
    bool HasChannel(const std::string &channel_name) const;

    bool HasChannel(int channel_number) const
    { return channel_number >= 0 && channel_number < n_channels(); }

    /// if HasChannelName(), return channel_layouts_[channel_name]->HasPyrLevel
    bool HasChannelPyrLevel(const std::string& channel_name, int pyr_level) const;

    bool HasChannelPyrLevel(int channel_number, int pyr_level) const;

    /// assert constructed_ = true. call ReadChannelLayout for instances managed by channel_layouts_
    void RefreshChannelLayouts();

    /// assert constructed_ = true. assert HasChannelName(channel_name)
    /// call ReadChannelLayout for channel_layout_.at(channel_name)
    void RefreshChannelLayout(const std::string& channel_name);

    std::string volume_root_dir() const
    { return volume_root_dir_; }

    std::vector<std::string> channel_dir_names() const
    { return channel_dir_names_; }

    std::vector<std::string> channel_names() const
    { return channel_names_; }

    /// given a channel name, return its channel number. channel names are
    /// alphabetically sorted (ascending). return -1 if channel_name is not
    /// found in channel_names_. note that channel_dir_names_ will not be searched.
    /// if channel_dir_names_ has unique entries, it will be identical to
    /// channel_names_.
    int channel_number(const std::string& channel_name) const;

    /// return channel_names_[channel_number] if channel_number is in valid range
    std::string channel_name(int channel_number, int time = 0) const;

    std::string channel_dir_name(int channel_number, int time = 0) const;

    const MChannelLayout& channel_layout(const std::string& channel_name) const;

    const MChannelLayout& channel_layout(int channel_number, int time = 0) const
    { return channel_layout(channel_name(channel_number, time)); }

    std::string channel_root_dir(const std::string& channel_name, int time = 0) const;

    /// full path to given channel number
    std::string channel_root_dir(int channel_number, int time = 0) const
    { return mcp3d::JoinPath(volume_root_dir_, channel_dir_name(channel_number, time)); }

    int n_channels() const
    { return (int)(channel_names_.size()); }

    int n_times() const
    { return 1; }

    friend class MImageInfo;

private:
    /// called by instances associated with non in memory images
    /// set channel_dir_names_ to channel_dir_names. if channel_dir_names_
    /// is empty, set channel_dir_names to {""}
    /// verify that all channel directories (join volume_root_ and
    /// channel_dir_names entries) exist. throw exception if non existent
    /// channel directory encountered
    /// set channel_names_ to channel_dir_names_
    void set_channel_names(const std::vector<std::string> &channel_dir_names);

    std::string channel_dir_name_strings() const;

    /// construct ChannelLayout for volume_root_ + channel_dir_names_
    /// if imaris file format found under channel layout, assert that a single
    /// empty string exists in channel_dir_names_. reset channel_dir_names,
    /// channel_names_ and channel_layouts_ according to channel numbers and
    /// channel names of the imaris file
    void InitChannelLayouts();

    /// return if empty()
    /// (1) channel_names_, channel_dir_names_, channel_layouts_ have the same size()
    /// (2) channel_names_ entries unique
    /// (3) all channel_names_ entry in channel_layouts_
    /// (4) if one channel has IMARIS format, all channels have IMARIS format
    /// (5) channel_layouts_.at(channel_name)->AssertValid()
    void AssertValid() const;

    const std::shared_ptr<MChannelLayout>& channel_layout_ptr(const std::string &channel_name) const;

    const std::string volume_root_dir_;
    // channel_dir_names_ is the name of the directory where data for the channel
    // is the directory of channel i data. channel_names_ is the name of each channel.
    // for formats other than imaris, channel_names_[i] is same as channel_dir_names_[i].
    // for imaris format, channel_dir_names_ are group names for the channel data
    // within the imaris file
    // is located. for formats other than imaris, volume_root_dir_/channel_dir_names_[i]
    std::vector<std::string> channel_dir_names_, channel_names_;
    // channel_name: MChannelLayout
    std::unordered_map<std::string, std::shared_ptr<MChannelLayout>> channel_layouts_;
    bool constructed_;
};

}

#endif //MCP3D_MCP3D_PROJECT_LAYOUT_HPP
