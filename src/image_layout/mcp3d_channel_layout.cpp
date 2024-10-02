//
// Created by muyezhu on 3/23/18.
//
using namespace std;

#include <regex>
#include <algorithm>
#include <boost/filesystem.hpp>
#include "image_interface/mcp3d_imaris_util.hpp"
#include "mcp3d_channel_layout.hpp"

mcp3d::MChannelLayout::MChannelLayout(const string& channel_root_dir): channel_root_dir_(channel_root_dir)
{
    if (!mcp3d::IsDir(channel_root_dir_))
        MCP3D_OS_ERROR(channel_root_dir_ + " is not a directory")
    ReadChannelLayout();
}

mcp3d::MChannelLayout::MChannelLayout(const mcp3d::MChannelLayout& other):
        channel_root_dir_(other.channel_root_dir_), channel_pyr_slices_(unordered_map<int, shared_ptr<mcp3d::MChannelPyrSlices>>{})
{
    for (const auto& item: other.channel_pyr_slices_)
        channel_pyr_slices_[item.first] = make_shared<mcp3d::MChannelPyrSlices>(*(item.second.get()));
}

void mcp3d::MChannelLayout::ReadChannelLayout()
{
    ClearChannelPyrSlices();
    // read/refresh channel pyr slices structure for pyr levels currently found under channel_root_dir_
    for (const auto& pyr_dir: PyrLevelDirs(true))
    {
        int pyr_level = mcp3d::MChannelLayout::DirPyrLevel(pyr_dir);
        if (channel_pyr_slices_.find(pyr_level) == channel_pyr_slices_.end())
            channel_pyr_slices_[pyr_level] = make_shared<mcp3d::MChannelPyrSlices>(pyr_dir);
        else
            // call MChannelPyrSlices::set_channel_pyr_dir since channel_pyr_dir_ has been cleared
            channel_pyr_slices_[pyr_level]->set_channel_pyr_dir(pyr_dir);
    }
    if (channel_pyr_slices_[0]->file_format() == mcp3d::FileFormat::IMARIS)
    {
        MCP3D_ASSERT(channel_pyr_slices_[0]->channel_pyr_dir() == channel_root_dir_)
        string imaris_path = mcp3d::StitchedImarisPathInDir(channel_pyr_slices_[0]->channel_pyr_dir());
        int n = mcp3d::NumberOfImarisResolutions(imaris_path);
        // if channel_pyr_slices_[i] does not previously exist, construct from channel_root_dir_,
        // otherwise call MChannelPyrSlices::set_channel_pyr_dir
        for (int i = 1; i < n; ++i)
        {
            if (channel_pyr_slices_.find(i) == channel_pyr_slices_.end())
                channel_pyr_slices_[i] = make_shared<mcp3d::MChannelPyrSlices>(channel_root_dir_);
            else
                channel_pyr_slices_[i]->set_channel_pyr_dir(channel_root_dir_);
        }
    }
    // for pyr levels no longer existing on disk, call MChannelPyrSlices::Clear
    for (auto& item: channel_pyr_slices_)
    {
        if (!mcp3d::IsDir(item.second->channel_pyr_dir()))
            item.second->Clear();
    }
}

bool mcp3d::MChannelLayout::operator==(const mcp3d::MChannelLayout& other) const
{
    if (channel_root_dir_ != other.channel_root_dir_)
        return false;
    if (channel_pyr_slices_.size() != other.channel_pyr_slices_.size())
        return false;
    for (const auto& item: channel_pyr_slices_)
    {
        if (other.channel_pyr_slices_.find(item.first) == other.channel_pyr_slices_.end())
        {
            if (!item.second->empty())
                return false;
        }
        else
        {
            if (*(item.second.get()) != *(other.channel_pyr_slices_.at(item.first).get()))
                return false;
        }
    }
    for (const auto& item: other.channel_pyr_slices_)
    {
        if (channel_pyr_slices_.find(item.first) == channel_pyr_slices_.end())
        {
            if (!item.second->empty())
                return false;
        }
    }
    return true;
}

bool mcp3d::MChannelLayout::HasPyrLevel(int pyr_level) const
{
    if (channel_pyr_slices_.find(pyr_level) == channel_pyr_slices_.end())
        return false;
    else
        return !(channel_pyr_slices_.at(pyr_level)->empty());
}

vector<string> mcp3d::MChannelLayout::PyrLevelDirs(bool full_path) const
{
    vector<string> pyr_level_dirs(1, full_path ? channel_root_dir_ : "");
    boost::filesystem::directory_iterator dir_iter{channel_root_dir_};
    for(auto it = boost::filesystem::begin(dir_iter); it != boost::filesystem::end(dir_iter); ++it)
    {
        if (boost::filesystem::is_directory(it->path()))
        {
            string dir_name = it->path().filename().string();
            if (mcp3d::MChannelLayout::IsPyrLevelDir(dir_name))
                pyr_level_dirs.push_back(full_path ? mcp3d::JoinPath({channel_root_dir_, dir_name}) : dir_name);
        }
    }
    sort(pyr_level_dirs.begin(), pyr_level_dirs.end());
    return pyr_level_dirs;
}


vector<int> mcp3d::MChannelLayout::pyr_levels() const
{
    vector<int> levels;
    for (const auto& item: channel_pyr_slices_)
    {
        if (!item.second->empty())
            levels.push_back(item.first);
    }
    sort(levels.begin(), levels.end());
    return levels;
}

string mcp3d::MChannelLayout::PyrLevelDirName(int pyr_level)
{
    MCP3D_ASSERT(pyr_level >= 0 && pyr_level < 99)
    if (pyr_level == 0)
        return string{};
    else
        return pyr_level_dir_prefix().append(mcp3d::PadNumStr(pyr_level, 2));
}

bool mcp3d::MChannelLayout::IsPyrLevelDir(const string &dir)
{
    string dir_name = mcp3d::Basename(dir);
    return regex_match(dir_name, mcp3d::MChannelLayout::pyr_level_dir_pattern());
}

int mcp3d::MChannelLayout::DirPyrLevel(const std::string &dir)
{
    string dir_name = mcp3d::Basename(dir);
    smatch m;
    if (!regex_match(dir_name, m, mcp3d::MChannelLayout::pyr_level_dir_pattern()))
        return 0;
    return stoi(m.str(1));
}

mcp3d::FileFormat mcp3d::MChannelLayout::file_format(int pyr_level) const
{
    MCP3D_ASSERT(HasPyrLevel(pyr_level))
    return channel_pyr_slices_.at(pyr_level)->file_format();
}

string mcp3d::MChannelLayout::channel_pyr_dir(int pyr_level)
{
    if (!HasPyrLevel(pyr_level))
        ReadChannelLayout();
    MCP3D_ASSERT(HasPyrLevel(pyr_level))
    return mcp3d::JoinPath(channel_root_dir_, mcp3d::MChannelLayout::PyrLevelDirName(pyr_level));
}

const shared_ptr<mcp3d::MChannelPyrSlices>& mcp3d::MChannelLayout::channel_pyr_slices_ptr(int pyr_level) const
{
    MCP3D_ASSERT(HasPyrLevel(pyr_level))
    return channel_pyr_slices_.at(pyr_level);
}

void mcp3d::MChannelLayout::AssertValid() const
{
    MCP3D_ASSERT(HasPyrLevel(0))
    for (const auto& item: channel_pyr_slices_)
        MCP3D_ASSERT(item.second)
    if (channel_pyr_slices_.at(0)->file_format() == mcp3d::FileFormat::IMARIS)
    {
        MCP3D_ASSERT(NumberOfPyrLevelDirs() == 1)
        string imaris_path = mcp3d::StitchedImarisPathInDir(channel_pyr_slices_.at(0)->channel_pyr_dir());
        int n = mcp3d::NumberOfImarisResolutions(imaris_path);
        vector<int> levels = pyr_levels();
        MCP3D_ASSERT(n == (int)levels.size())
        for (int i = 1; i < n; ++i)
        {
            MCP3D_ASSERT(HasPyrLevel(i))
            MCP3D_ASSERT(channel_pyr_slices_.at(i)->channel_pyr_dir() == channel_pyr_slices_.at(0)->channel_pyr_dir())
            MCP3D_ASSERT(channel_pyr_slices_.at(i)->file_format() == mcp3d::FileFormat::IMARIS)
        }
    }
    else
    {
        for (const auto& item: channel_pyr_slices_)
            MCP3D_ASSERT(item.second->file_format() != mcp3d::FileFormat::IMARIS)
        for (const auto& item: channel_pyr_slices_)
            MCP3D_ASSERT(HasPyrLevel(item.first) == HasPyrLevelDir(item.first))
        for (const auto& pyr_level_dir: PyrLevelDirs())
        {
            int pyr_level = MChannelLayout::DirPyrLevel(pyr_level_dir);
            MCP3D_ASSERT(HasPyrLevel(pyr_level) && HasPyrLevelDir(pyr_level))
        }
    }
}

void mcp3d::MChannelLayout::ClearChannelPyrSlices()
{
    for (auto& item: channel_pyr_slices_)
        item.second->Clear();
}

