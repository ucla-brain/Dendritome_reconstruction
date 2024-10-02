//
// Created by muyezhu on 3/2/19.
//
#include <string>
#include <regex>
#include <utility>
#include <unordered_set>
#include <algorithm>
#include "common/mcp3d_common.hpp"
#include "image_interface/mcp3d_imaris_util.hpp"
#include "mcp3d_image_layout_constants.hpp"
#include "mcp3d_volume_layout.hpp"

using namespace std;

mcp3d::MVolumeLayout::MVolumeLayout(const string &volume_root_dir, const vector<string>& channel_dir_names):
                           volume_root_dir_(volume_root_dir), constructed_(false)
{
    if (!mcp3d::IsDir(volume_root_dir_))
       MCP3D_RUNTIME_ERROR(volume_root_dir_ + " is not a directory")
    set_channel_names(channel_dir_names);
    InitChannelLayouts();
    constructed_ = true;
}

mcp3d::MVolumeLayout::MVolumeLayout(const MVolumeLayout& other):
        volume_root_dir_(other.volume_root_dir_), channel_dir_names_(other.channel_dir_names_), channel_names_(other.channel_names_), constructed_(other.constructed_)
{
    for (const auto& item: other.channel_layouts_)
        channel_layouts_[item.first] = make_shared<mcp3d::MChannelLayout>(*(item.second.get()));
    MCP3D_ASSERT(constructed_)
}

void mcp3d::MVolumeLayout::Print() const
{
    cout << "volume root = " << volume_root_dir_ << endl;
    cout << "channel directory names = " << channel_dir_name_strings() << endl;
}

bool mcp3d::MVolumeLayout::operator==(const MVolumeLayout &other) const
{
    if (volume_root_dir_ != other.volume_root_dir_)
        return false;
    if (channel_dir_names_ != other.channel_dir_names_ || channel_names_ != other.channel_names_)
        return false;
    if (channel_layouts_.size() != other.channel_layouts_.size())
        return false;
    for (const auto& item: channel_layouts_)
    {
        if (*(item.second.get()) != *(other.channel_layouts_.at(item.first).get()))
            return false;
    }
    return true;
}

bool mcp3d::MVolumeLayout::HasChannel(const string &channel_name) const
{
    for (const auto& name: channel_names_)
        if (name == channel_name)
            return true;
    return false;
}

bool mcp3d::MVolumeLayout::HasChannelPyrLevel(const std::string &channel_name, int pyr_level) const
{
    if (!HasChannel(channel_name))
        return false;
    return channel_layouts_.at(channel_name)->HasPyrLevel(pyr_level);
}

bool mcp3d::MVolumeLayout::HasChannelPyrLevel(int channel_number, int pyr_level) const
{
    if (!HasChannel(channel_number))
        return false;
    return channel_layouts_.at(channel_name(channel_number, 0))->HasPyrLevel(pyr_level);
}

void mcp3d::MVolumeLayout::RefreshChannelLayouts()
{
    MCP3D_ASSERT(constructed_)
    for (auto& item: channel_layouts_)
        item.second->ReadChannelLayout();
}

void mcp3d::MVolumeLayout::RefreshChannelLayout(const std::string &channel_name)
{
    MCP3D_ASSERT(constructed_)
    if (!HasChannel(channel_name))
        MCP3D_RUNTIME_ERROR("channel name " + channel_name + " is not found under the MVolumeLayout instance (volume_root_dir = " + volume_root_dir_ + ")")
    channel_layouts_.at(channel_name)->ReadChannelLayout();
}

string mcp3d::MVolumeLayout::channel_name(int channel_number, int time) const
{
    if (!HasChannel(channel_number))
        MCP3D_RUNTIME_ERROR("channel number " + to_string(channel_number) + " is not found under the MVolumeLayout instance (volume_root_dir = " + volume_root_dir_ + ")")
    MCP3D_ASSERT(time == 0)
    return channel_names_[channel_number];
}

string mcp3d::MVolumeLayout::channel_dir_name(int channel_number, int time) const
{
    if (!HasChannel(channel_number))
        MCP3D_RUNTIME_ERROR("channel number " + to_string(channel_number) + " is not found under the MVolumeLayout instance (volume_root_dir = " + volume_root_dir_ + ")")
    MCP3D_ASSERT(time == 0)
    return channel_dir_names_[channel_number];
}

string mcp3d::MVolumeLayout::channel_dir_name_strings() const
{
    string output;
    for (const auto& channel_dir_name: channel_dir_names_)
    {
        if (!channel_dir_name.empty())
            output.append(channel_dir_name);
        else
            output.append("\"\"");
        output.append("    ");
    }
    return output;
}

const mcp3d::MChannelLayout& mcp3d::MVolumeLayout::channel_layout(const std::string &channel_name) const
{
    if (!HasChannel(channel_name))
        MCP3D_RUNTIME_ERROR("channel name " + channel_name + " is not found under the MVolumeLayout instance (volume_root_dir = " + volume_root_dir_ + ")")
    return *(channel_layouts_.at(channel_name).get());
}

string mcp3d::MVolumeLayout::channel_root_dir(const std::string &channel_name, int time) const
{
    if (!HasChannel(channel_name))
        MCP3D_RUNTIME_ERROR("channel name " + channel_name + " is not found under the MVolumeLayout instance (volume_root_dir = " + volume_root_dir_ + ")")
    MCP3D_ASSERT(time == 0)
    return mcp3d::JoinPath(volume_root_dir_, channel_dir_name(channel_number(channel_name)));
}

const shared_ptr<mcp3d::MChannelLayout>& mcp3d::MVolumeLayout::channel_layout_ptr(const string &channel_name) const
{
    if (channel_layouts_.find(channel_name) == channel_layouts_.end())
        MCP3D_RUNTIME_ERROR("channel name " + channel_name + " is not found under the MVolumeLayout instance (volume_root_dir = " + volume_root_dir_ + ")")
    return channel_layouts_.at(channel_name);
}

int mcp3d::MVolumeLayout::channel_number(const string &channel_name) const
{
    for (int i = 0; i < n_channels(); ++i)
    {
        if (channel_names_[i] == channel_name)
            return i;
    }
    return -1;
}

void mcp3d::MVolumeLayout::set_channel_names(const vector<string> &channel_dir_names)
{
    MCP3D_ASSERT(!constructed_)
    channel_dir_names_ = channel_dir_names;
    if (channel_dir_names_.empty())
        channel_dir_names_ = vector<string>({""});
    sort(channel_dir_names_.begin(), channel_dir_names_.end());
    for (const auto& channel_dir_name: channel_dir_names_)
    {
        string channel_dir = mcp3d::JoinPath(volume_root_dir_, channel_dir_name);
        if (!mcp3d::IsDir(channel_dir))
            MCP3D_RUNTIME_ERROR(channel_dir + " is not a directory")
    }
    channel_names_ = channel_dir_names_;
}

void mcp3d::MVolumeLayout::InitChannelLayouts()
{
    MCP3D_ASSERT(!constructed_)
    for (int i = 0; i < n_channels(); ++i)
    {
        channel_layouts_[channel_names_[i]] = make_shared<mcp3d::MChannelLayout>(channel_root_dir(i, 0));
        if (channel_layouts_[channel_names_[i]]->file_format(0) == mcp3d::FileFormat::IMARIS)
        {
            // if imaris file is found, there should be only one channel folder directory
            MCP3D_ASSERT(n_channels() == 1)
            // if imaris file is found, channel dir name should be empty string
            MCP3D_ASSERT(channel_dir_names_[0].empty())
            // fill channel_names_, channel_dir_names_ and channel_layouts_ with the imaris file
            string imaris_path = mcp3d::StitchedImarisPathInDir(volume_root_dir_);
            channel_names_ = mcp3d::ImarisResolutionChannelNames(imaris_path, 0);
            channel_dir_names_ = vector<string>((int)channel_names_.size(), "");
            channel_layouts_.clear();
            for (const auto& channel_name_: channel_names_)
                channel_layouts_[channel_name_] = make_shared<mcp3d::MChannelLayout>(volume_root_dir_);
        }
    }
}

void mcp3d::MVolumeLayout::AssertValid() const
{
    if (empty())
        return;
    MCP3D_ASSERT(channel_names_.size() == channel_dir_names_.size())
    MCP3D_ASSERT(channel_names_.size() == channel_layouts_.size())
    MCP3D_ASSERT(channel_names_.size() == unordered_set<string>(channel_names_.begin(), channel_names_.end()).size())
    bool is_imaris = false;
    for (const auto& channel_name: channel_names_)
    {
        MCP3D_ASSERT(channel_layouts_.find(channel_name) != channel_layouts_.end())
        channel_layouts_.at(channel_name)->AssertValid();
        if (channel_layouts_.at(channel_name)->file_format(0) == mcp3d::FileFormat::IMARIS)
            is_imaris = true;
    }
    if (is_imaris)
    {
        for (const auto& channel_name: channel_names_)
            MCP3D_ASSERT(channel_layouts_.at(channel_name)->file_format(0) == mcp3d::FileFormat::IMARIS)
    }
}
