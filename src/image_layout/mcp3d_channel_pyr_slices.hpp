//
// Created by muyezhu on 3/17/19.
//

#ifndef MCP3D_MCP3D_FOLDER_SLICES_HPP
#define MCP3D_MCP3D_FOLDER_SLICES_HPP

#include <memory>
#include <regex>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <nlohmann/json.hpp>
#include "common/mcp3d_common.hpp"
#include "image_interface/mcp3d_file_formats.hpp"
#include "mcp3d_image_layout_constants.hpp"

/* an image having multiple channels have number of xyz volumes equal to number
 * of channels. under each channel, folder slice layout, pyr level layout are
 * allowed to be different
 * */
namespace mcp3d
{
/// one of the 3 dimensions that a given folder slice into
enum class ChannelAxis
{
    X = 0, Y = 1, Z = 2
};

inline std::string ChannelAxisStr(ChannelAxis axis)
{
    if (axis == ChannelAxis::X)
        return "x";
    else if (axis == ChannelAxis::Y)
        return "y";
    else
        return "z";
}

inline ChannelAxis StrToChannelAxis(const std::string& axis_str)
{
    if (axis_str == "z")
        return ChannelAxis::Z;
    else if (axis_str == "y")
        return ChannelAxis::Y;
    else if (axis_str == "x")
        return ChannelAxis::X;
    else
        MCP3D_INVALID_ARGUMENT("axis_str must be one of \"z\", \"y\" or \"x\"");
}

/// MAX_AXIS_SLICE_DIM: represents the slice dimension when there is a single slice along an axis
/// the corresponding slice_numbers_ entry along the axis should be 1
class MChannelPyrSlices
{
public:
    /// default constructor creates an empty instance
    MChannelPyrSlices(): slice_dims_(std::map<ChannelAxis, int>{}), slice_numbers_(std::map<ChannelAxis, int>{}),
                         channel_pyr_dir_(std::string{}), file_format_(FileFormat::UNKNOWN) {}

    explicit MChannelPyrSlices(const std::string& channel_pyr_dir)
    { set_channel_pyr_dir(channel_pyr_dir); }

    MChannelPyrSlices(const MChannelPyrSlices& other) = default;

    bool operator==(const MChannelPyrSlices& other) const;

    bool operator!=(const MChannelPyrSlices& other) const
    { return !(*this == other); }

    bool empty() const
    { return channel_pyr_dir_.empty(); }

    /// if channel_pyr_dir_ does not exist (any more), clear instance and return
    /// otherwise read folders and file formats
    void Refresh();

    /// set instance to empty state
    void Clear();

    /// create folder structures according to parameters.
    /// channel_pyr_dir will be made or deleted and remade as an empty directory. sub directories are then populated
    /// file_format_ is unknown after Make, since all leaf sub directories will be empty
    void Make(const std::string& channel_pyr_dir, const std::map<ChannelAxis, int> slice_numbers = std::map<ChannelAxis, int>{},
              const std::map<ChannelAxis, int> slice_dims = std::map<ChannelAxis, int>{});

    /// regex pattern of valid slice name e.g. "z0000_001024/y0001_000256/x0001_000256"
    /// this does not cover "", which validly represents cases where there's a
    /// single slice covering the entirety of an axis
    static std::regex AxisSliceNamePattern()
    { return std::regex{"^([xyz])([0-9]{4})_([0-9]{6})$"}; }

    /// if slice_dim equal to MAX_AXIS_SLICE_DIM, return empty string
    /// assert slice_dim and slice_id to be non negative and within maximum width
    static std::string AxisSliceName(ChannelAxis axis, int slice_id, int slice_dim);

    /// return full slice name specified by the slice ids by joining outputs from AxisSliceNameFromSliceIDs
    /// e.g. "z0000_001024/y0001_000256/x0001_000256"
    std::string SliceNameFromSliceIDs(int zslice_id, int yslice_id, int xslice_id) const;

    /// return full slice name that will contain the voxel (z_coordinate, y_coordinate, x_coordinate), by joining
    /// outputs from AxisSliceNameFromCoordinates
    /// the coordinates here are local to the current pyr_level volume
    std::string SliceNameFromCoordinates(int z_coordinate, int y_coordinate, int x_coordinate) const;

    std::vector<std::string> SliceNames() const;

    /// return false if empty
    /// return true if the channel pyramid folder is flat (single slice along every axis)
    bool IsFlat() const
    { return AxisIsFlat(ChannelAxis::Z) && AxisIsFlat(ChannelAxis::Y) && AxisIsFlat(ChannelAxis::X); }

    /// return false if empty
    bool AxisIsFlat(ChannelAxis axis) const
    { return empty() ? false : slice_dims_.at(axis) == MAX_AXIS_SLICE_DIM; }

    /// return -1 if empty (0 is MAX_AXIS_SLICE_DIM constant)
    int slice_dim(ChannelAxis axis) const
    { return empty() ? -1 : slice_dims_.at(axis); }

    /// return 0 if empty
    int slice_number(ChannelAxis axis) const
    { return empty() ? 0 : slice_numbers_.at(axis); }

    std::string channel_pyr_dir() const
    { return channel_pyr_dir_; }

    FileFormat file_format() const
    { return file_format_; }

    /// assert channel_pyr_dir is directory. clear instance. set channel_pyr_dir_
    /// to channel_pyr_dir. call Refresh
    void set_channel_pyr_dir(const std::string& channel_pyr_dir);

private:
    /// if empty() return empty string
    /// "[axis][slice_id]_[slice_dim]" for xyz slices, e.g. x0001_001024.
    /// for slices that cover the entire axis, return empty string
    /// assert slice_id within valid range
    std::string AxisSliceNameFromSliceID(ChannelAxis axis, int slice_id) const;

    /// if empty() return empty string
    /// for slices that cover the entire axis, return empty string
    /// otherwise assert coordinate within valid range
    /// return "[axis][slice_id]_[slice_dim]" for xyz slices that contains the given coordinate
    /// e.g. if along x axis, slice_numbers = 2, slice_dim = 1024, then coordinate 2000
    /// is contained in x0001_1024.
    std::string AxisSliceNameFromCoordinate(ChannelAxis axis, int coordinate) const;

    /// if empty instance, do nothing. otherwise assert channel_pyr_dir_ exists
    /// clears slice_dim_, slice_numbers_. set file_format_ to unknown. read folder slices and contained file format under channel_pyr_dir_
    void ReadChannelPyrSlices();

    void ReadChannelPyrSlicesTree(const std::string& input_dir, std::unordered_map<std::string, std::unordered_set<int>>& slice_ids);

    /// if empty instance, do nothing. otherwise assert channel_pyr_dir_ exists.
    /// search order: FileFormatSearchOrder().
    /// for each format, first determine if channel_pyr_dir_'s directory layout support the format.
    /// if true, scan all slice directories for files that are of the format
    /// set file_format_ to the first format found. if not such format exists, file_format_ is UNKNOWN
    void ReadChannelPyrFormat();

    /// if format is UNKNOWN, return false. empty() instances should return false for any format
    bool ChannelPyrDirCompatibleWithFileFormat(FileFormat file_format);

    /// return true if folder slices are flat along xy axes
    bool ChannelPyrDirCompatibleWithTiffFormat() const
    { return AxisIsFlat(ChannelAxis::Y) && AxisIsFlat(ChannelAxis::X); }

    bool ChannelPyrDirCompatibleWithOmeTiffFormat() const
    { return !empty(); }

    /// return true if channel_pyr_dir_ is flat along all zyx axes.
    /// MChannelLayout will further assert that no other pyr level dirs are present in channel_pyr_dir.
    bool ChannelPyrDirCompatibleWithImarisFormat() const
    { return IsFlat(); }

    bool ChannelPyrDirCompatibleWithHdf5Format() const
    { return !empty(); }

    /// slice dimensions is local to the pyr_level. it has no awareness of
    /// the pyramid mapping or the global image volume dimensions
    std::map<ChannelAxis, int> slice_dims_, slice_numbers_;
    std::string channel_pyr_dir_;
    FileFormat file_format_;
};

}

#endif //MCP3D_MCP3D_FOLDER_SLICES_HPP
