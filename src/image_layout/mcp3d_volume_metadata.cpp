//
// Created by muyezhu on 8/19/20.
//

#include <regex>
#include <vector>
#include "common/mcp3d_common.hpp"
#include "mcp3d_volume_metadata.hpp"

using namespace std;

mcp3d::MVolumeMetadata::MVolumeMetadata(const string &volume_root_dir):
        volume_root_dir_(volume_root_dir), scope_type_(ScopeType::UNKNOWN), tissue_type_(TissueType::UNKNOWN),
        cutting_plane_(CuttingPlane::UNKNOWN), magnification_(0), thickness_(0), slice_name_(std::string{}),
        slide_name_(std::string{}), region_name_(std::string{})
{
    if (!mcp3d::IsDir(volume_root_dir_))
        MCP3D_OS_ERROR(volume_root_dir_ + " is not a directory")
    ParseMetaData();
}

void mcp3d::MVolumeMetadata::ParseMetaData()
{
    vector<string> path_components = mcp3d::PathComponents(volume_root_dir_);
    size_t component_id = 0;
    // get scope type and magnification from path components
    // if found, continue to tissue type
    while (component_id < path_components.size())
    {
        pair<ScopeType, int> scope_mag = mcp3d::MVolumeMetadata::ParseScopeMagnification(path_components[component_id]);
        ++component_id;
        if (scope_mag.first != mcp3d::ScopeType::UNKNOWN)
        {
            scope_type_ = scope_mag.first;
            magnification_ = scope_mag.second;
            break;
        }
    }

    if (component_id < path_components.size())
    {
        // parse tissue type
        tissue_type_ = mcp3d::MVolumeMetadata::ParseTissueType(path_components[component_id]);
        ++component_id;
        // parse slice
        if (component_id < path_components.size())
            ParseSlice(path_components[component_id]);
    }
}

void mcp3d::MVolumeMetadata::ParseSlice(const std::string &path_component)
{
    // parse slice+
    // A repeated capturing group will only capture the last iteration.
    // Put a capturing group around the repeated group to capture all iterations
    // or use a non-capturing group instead if you're not interested in the data
    regex slice_pattern("^([0-9]+_[0-9]+)(_[0-9]+um)*((_[a-zA-Z]+)*)$");
    smatch match;
    bool has_slice = regex_match(path_component, match, slice_pattern);
    if (has_slice)
    {
        slice_name_ = match.str(0);
        slide_name_ = match.str(1);
        if (!match.str(2).empty())
        {
            string thickness_str = match.str(2).substr(1, match.str(2).size() - 3);
            thickness_ = stoi(thickness_str);
        }
        if (!match.str(3).empty())
            region_name_ = match.str(3).substr(1);
        if (!match.str(4).empty())
        {
            region_name_ = match.str(3).substr(1, match.str(3).size() - match.str(4).size() - 1);
            cutting_plane_ = mcp3d::MVolumeMetadata::ParseCuttingPlane(match.str(4).substr(1));
        }
    }
}

pair<mcp3d::ScopeType, int> mcp3d::MVolumeMetadata::ParseScopeMagnification(const string& path_component)
{
    ScopeType scope_type;
    int mag;
    for (const auto& scope_mag: mcp3d::ValidScopeMagnifications)
    {
        scope_type = scope_mag.first;
        mag = scope_mag.second;
        string result(mcp3d::MVolumeMetadata::ScopeMagnificationString(scope_type, mag));
        if (mcp3d::StringLower(result) == mcp3d::StringLower(path_component))
            return make_pair(scope_type, mag);
    }
    return make_pair(mcp3d::ScopeType::UNKNOWN, 0);
};

mcp3d::TissueType mcp3d::MVolumeMetadata::ParseTissueType(const std::string &path_component)
{
    string input(mcp3d::StringLower(path_component));
    if (input == "wholebrain" || input == "whole brain")
        return mcp3d::TissueType::WHOLE_BRAIN;
    else if (input == "lefthem" || input == "lefthemisphere" ||
             input == "left hemisphere" || input == "left hem")
        return mcp3d::TissueType::LEFT_HEMISPHERE;
    else if (input == "righthem" || input == "righthemisphere" ||
             input == "right hemisphere" || input == "right hem")
        return mcp3d::TissueType::RIGHT_HEMISPHERE;
    else if (input == "spinalcord" || input == "spinal cord")
        return mcp3d::TissueType::SPINAL_CORD;
    else
        return mcp3d::TissueType::UNKNOWN;
}

mcp3d::CuttingPlane mcp3d::MVolumeMetadata::ParseCuttingPlane(const std::string &cutting_plane_string)
{
    if (mcp3d::StringLower(cutting_plane_string) == "coronal")
        return mcp3d::CuttingPlane::CORONAL;
    else if (mcp3d::StringLower(cutting_plane_string) == "saggital")
        return mcp3d::CuttingPlane::SAGGITAL;
    else if (mcp3d::StringLower(cutting_plane_string) == "none")
        return mcp3d::CuttingPlane::NONE;
    else
        return mcp3d::CuttingPlane::UNKNOWN;
}

string mcp3d::MVolumeMetadata::ScopeTypeString(mcp3d::ScopeType scope_type)
{
    if (scope_type == ScopeType::DRAGONFLY)
        return "Dragonfly";
    else if (scope_type == ScopeType::LIGHTSHEET)
        return "LightSheet";
    else if (scope_type == ScopeType::FMOST)
        return "fMost";
    else
        return "unknown";
}

string mcp3d::MVolumeMetadata::ScopeMagnificationString(mcp3d::ScopeType scope_type, int mag)
{
    string result(mcp3d::MVolumeMetadata::ScopeTypeString(scope_type));
    result.append("_");
    result.append(to_string(mag));
    result.append("x");
    return result;
}


