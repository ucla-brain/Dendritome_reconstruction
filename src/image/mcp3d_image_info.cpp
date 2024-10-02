//
// Created by muyezhu on 6/20/19.
//
#include <utility>
#include <tuple>
#include <algorithm>
#include "common/mcp3d_utility.hpp"
#include "image_layout/mcp3d_volume_layout.hpp"
#include "mcp3d_image_constants.hpp"
#include "mcp3d_image_io.hpp"

using namespace std;

mcp3d::MImageInfo::MImageInfo(const string &volume_path, const vector<std::string> &channel_dir_names):
        volume_layout_(volume_path, channel_dir_names), channel_infos_(std::unordered_map<std::string, MChannelInfo>{})
{
    InitChannelInfos();
    AssertValid();
}

void mcp3d::MImageInfo::InitChannelInfos()
{
    for (const auto& channel_name: volume_layout_.channel_names())
    {
        mcp3d::MChannelInfo channel_info(volume_layout_.channel_layout_ptr(channel_name));
        channel_infos_.emplace(channel_name, move(channel_info));
    }
}

void mcp3d::MImageInfo::AddChannelPyrInfo(const string& channel_name, mcp3d::MChannelPyrInfo&& channel_pyr_info)
{
    if (!volume_layout_.HasChannel(channel_name))
        MCP3D_RUNTIME_ERROR("channel name " + channel_name + " not found under " + volume_layout_.volume_root_dir())
    channel_infos_.at(channel_name).AddPyrInfo(move(channel_pyr_info));
    AssertValid();
}

void mcp3d::MImageInfo::Save()
{
    AssertValid();
    for (const auto& item: channel_infos_)
    {
        item.second.Save();
        if (item.second.file_format(0) == mcp3d::FileFormat::IMARIS)
            break;
    }
}

bool mcp3d::MImageInfo::HasChannelPyrInfo(const string &channel_name, int pyr_level) const
{
    if (!volume_layout_.HasChannel(channel_name))
        return false;
    return channel_infos_.at(channel_name).HasPyrInfo(pyr_level);
}

bool mcp3d::MImageInfo::HasChannelPyrInfo(int channel_number, int pyr_level) const
{
    if (!volume_layout_.HasChannel(channel_number))
        return false;
    return channel_infos_.at(volume_layout_.channel_name(channel_number, 0)).HasPyrInfo(pyr_level);
}

int mcp3d::MImageInfo::n_pyr_levels(const string &channel_name) const
{
    if (!volume_layout_.HasChannel(channel_name))
        MCP3D_RUNTIME_ERROR(channel_name + " not found under " + volume_layout_.volume_root_dir())
    return channel_infos_.at(channel_name).n_pyr_levels();
}

int mcp3d::MImageInfo::n_pyr_infos(const string &channel_name) const
{
    if (!volume_layout_.HasChannel(channel_name))
        MCP3D_RUNTIME_ERROR(channel_name + " not found under " + volume_layout_.volume_root_dir())
    return channel_infos_.at(channel_name).n_pyr_infos();
}

int mcp3d::MImageInfo::MaxMappablePyrLevel() const
{
    for (const auto& channel_name: volume_layout_.channel_names())
    {
        if (channel_infos_.at(channel_name).n_pyr_infos() > 0)
        {
            for (const auto& item: channel_infos_.at(channel_name).channel_pyr_infos_)
            {
                int max_level = item.first, pyr_ratio = 1;
                int current_xdim = item.second.xdim(), current_ydim = item.second.ydim();
                while (true)
                {
                    if (current_xdim < pyr_ratio || current_ydim < pyr_ratio)
                        break;
                    pyr_ratio *= 2;
                    current_xdim /= 2;
                    current_ydim /= 2;
                    ++max_level;
                }
                return max_level - 1;
            }
        }
    }
    return -1;
}

vector<int> mcp3d::MImageInfo::xyz_dims(int pyr_level) const
{
    for (const auto& item: channel_infos_)
    {
        if (item.second.HasPyrLevel(pyr_level))
            return item.second.xyz_dims(pyr_level);
    }
    return vector<int>({0, 0, 0});
}

const mcp3d::MChannelInfo& mcp3d::MImageInfo::channel_info(const string& channel_name) const
{
    if (!volume_layout_.HasChannel(channel_name))
        MCP3D_RUNTIME_ERROR(channel_name + " not found under " + volume_layout_.volume_root_dir())
    return channel_infos_.at(channel_name);
}

mcp3d::MChannelInfo& mcp3d::MImageInfo::channel_info(const std::string &channel_name)
{
    if (!volume_layout_.HasChannel(channel_name))
        MCP3D_RUNTIME_ERROR(channel_name + " not found under " + volume_layout_.volume_root_dir())
    return channel_infos_.at(channel_name);
}

void mcp3d::MImageInfo::AssertValid() const
{
    if (empty())
        return;
    // assert volume layout valid
    volume_layout_.AssertValid();
    // channel_infos_ has and only has channel_names_ of volume_out_ as keys
    MCP3D_ASSERT(volume_layout_.n_channels() == (int)channel_infos_.size())
    for (const auto& channel_name: volume_layout_.channel_names())
    {
        MCP3D_ASSERT(channel_infos_.find(channel_name) != channel_infos_.end())
        // volume_layout_.channel_layouts_ store the same pointer as channel_infos_.channel_layout_
        MCP3D_ASSERT(volume_layout_.channel_layout_ptr(channel_name).get() == channel_infos_.at(channel_name).channel_layout_ptr().get())
        // assert channel_infos_ valid
        channel_infos_.at(channel_name).AssertValid();
    }
    // for MChannelInfo instances that all have MChannelPyrInfo at a given pyr_level,
    // assert the xyz dimensions of the level are identical
    int max_n_pyr_infos = 0;
    for (const auto& item: channel_infos_)
        max_n_pyr_infos = max(max_n_pyr_infos, item.second.pyr_info_levels()[item.second.n_pyr_infos() - 1]);
    for (int pyr_level = 0; pyr_level < max_n_pyr_infos; ++pyr_level)
    {
        vector<int> pyr_xyz_dims;
        for (const auto& item: channel_infos_)
        {
            if (item.second.HasPyrInfo(pyr_level))
            {
                if (pyr_xyz_dims.empty())
                    pyr_xyz_dims = item.second.xyz_dims(pyr_level);
                else
                    MCP3D_ASSERT(pyr_xyz_dims == item.second.xyz_dims(pyr_level));
            }
        }
    }
}



