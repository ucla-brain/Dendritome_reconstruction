//
// Created by muyezhu on 8/19/20.
//

#ifndef MCP3D_MCP3D_VOLUME_METADATA_HPP
#define MCP3D_MCP3D_VOLUME_METADATA_HPP

#include <string>
#include <utility>
#include <set>

namespace mcp3d
{

enum class CuttingPlane
{
    CORONAL, SAGGITAL, NONE, UNKNOWN = -1
};

enum class ScopeType
{
    DRAGONFLY, LIGHTSHEET, FMOST, UNKNOWN = -1
};

enum class TissueType
{
    LEFT_HEMISPHERE, RIGHT_HEMISPHERE, WHOLE_BRAIN, SPINAL_CORD, UNKNOWN = -1
};

static const std::set<std::pair<ScopeType, int>> ValidScopeMagnifications
        {std::make_pair(ScopeType::DRAGONFLY, 2),
         std::make_pair(ScopeType::DRAGONFLY, 10),
         std::make_pair(ScopeType::DRAGONFLY, 30),
         std::make_pair(ScopeType::LIGHTSHEET, 4),
         std::make_pair(ScopeType::LIGHTSHEET, 10),
         std::make_pair(ScopeType::LIGHTSHEET, 25),
         std::make_pair(ScopeType::FMOST, 40),
        };

/* Raw/Dragonfly_10x/WholeBrain/1_05_400um_BLA(_coronal)/(../volume_name/..)/.ims
 *     volume root = Raw/Dragonfly_10x/WholeBrain/1_05_400um_BLA(_coronal)/(../volume_name/..)
 * Raw/Lightsheet_10x/SpinalCord/(../(slice_name)/(volume_name)/..)/ex0_em0/.raw
 *     volume root = Raw/Lightsheet_10x/SpinalCord/(../(slice_name)/(volume_name)/..)
 * Raw/Fmost_40x/WholeBrain/CH1.tif
 *     volume root = Raw/Fmost_40x/WholeBrain
 */

/// the constructor will try to parse scope_type_, tissue_type_, slice_name_ etc
/// as well

/// scope type is not used to enforce volume layout structures. e.g., dragonfly
/// scope type directory can have multiple channel directories (if no stitched
/// imaris format is found). during pipeline processing, it's very likely
/// to produce non imaris file format outputs under a dragonfly scope type directory
class MVolumeMetadata
{
public:
    MVolumeMetadata(): volume_root_dir_(std::string{}), scope_type_(ScopeType::UNKNOWN), tissue_type_(TissueType::UNKNOWN),
                       cutting_plane_(CuttingPlane::UNKNOWN), magnification_(0), thickness_(0), slice_name_(std::string{}),
                       slide_name_(std::string{}), region_name_(std::string{}) {}

    MVolumeMetadata(const std::string& volume_root_dir);

    ScopeType scope_type() const  { return scope_type_; }

    TissueType tissue_type() const  { return tissue_type_; }

    CuttingPlane cutting_plane() const  { return cutting_plane_; }

    int magnification() const  { return magnification_; }

    int thickness() const  { return thickness_; }

    std::string slice_name() const  { return slice_name_; }

    std::string slide_name() const  { return slide_name_; }

    std::string region_name() const  { return region_name_; }

    /// if path_component equal to string representation of any entry in
    /// ValidScopeMagnifications, return corresponding pair<ScopeType, int>.
    /// otherwise return (ScopeType::UNKNOWN, 0)
    static std::pair<ScopeType, int> ParseScopeMagnification(const std::string &path_component);

    static CuttingPlane ParseCuttingPlane(const std::string& cutting_plane_string);

    static TissueType ParseTissueType(const std::string &path_component);

private:
    /// parse meta data member variables including scope type, magnification,
    /// tissue type, cutting plane, thickness, slice name, slide name, region name
    /// from volume_root_dir_
    void ParseMetaData();

    void ParseSlice(const std::string &path_component);

    static std::string ScopeTypeString(ScopeType scope_type);

    static std::string ScopeMagnificationString(ScopeType scope_type, int mag);

    const std::string volume_root_dir_;
    ScopeType scope_type_;
    TissueType tissue_type_;
    CuttingPlane cutting_plane_;
    int magnification_, thickness_;
    std::string slice_name_, slide_name_, region_name_;

};

}

#endif //MCP3D_MCP3D_VOLUME_METADATA_HPP
