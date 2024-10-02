//
// Created by muyezhu on 2/11/18.
//

#ifndef MCP3D_MCP3D_CHANNEL_INFO_HPP
#define MCP3D_MCP3D_CHANNEL_INFO_HPP

#include <memory>
#include <unordered_set>
#include <unordered_map>
#include <nlohmann/json.hpp>
#include "common/mcp3d_common.hpp"
#include "mcp3d_image_constants.hpp"
#include "image_interface/mcp3d_voxel_types.hpp"
#include "image_interface/mcp3d_file_formats.hpp"
#include "image_layout/mcp3d_channel_pyr_slices.hpp"
#include "image_layout/mcp3d_channel_layout.hpp"

/* a channel directory contains one channel of the image volume and the
 * downsampled pyramids of the channel, organized as following:
 * channel_root_dir:
 *    level 0 image file(s)
 *    (pyr_level_01)
 *    (pyr_level_02)
 *    ...
 *    (pyr_level_n(+))
 * MChannelInfo describe all data under channel_root_dir
 * MChannelPyrInfo describe data at a given pyramid level
 *
 * MChannelPyrInfo required json entries:
 *    format: FileFormat enum string
 *    x dimension: total size of image along x dimension across the entire image
 *                 volume at current pyramid level
 *    xchunk dimension: a chunk is the unit of storage for image data.
 *                      ome-tif and hdf5 may save volume tiles.
 *                      tif image sequence's unit of storage is a single tiff
 *                      image plane
 *                      within this storage unit, the size along x dimension is
 *                      xchunk dimension
 *    resolution level: the resolution level of the pyramid
 *    channel pyramid directory: full path of directory where current image pyramid
 *                               files are located.
 *    image sequences: base name of image chunks sorted by xyzct order.
 *                     to recove full path of an image in standalone MChannelPyrInfo:
 *                     image pyramid directory + image basename
 *    voxel type: VoxelType enum string
 *    pyramid ratio: factor of downsample from level 0 image. between level i
 *                   and level i + 1, the downsample factor is strictly 2. other
 *                   factors are invalid
 *                   f_i = x_i / x_0 = y_i / y_0
 *
 * MChannelInfo: collection of MChannelPyrInfo, sorted by ascending pyramid level
 *             find the nearest common directory from image pyramid directories,
 *             assign it as image root dir.
 */
namespace mcp3d
{

enum class VolumeComplete
{
    COMPLETE = 0, INCOMPLETE = 1, UNKNOWN = -1
};

/// xyz volume at a given pyramid level
class MChannelPyrInfo
{
public:
    MChannelPyrInfo(): channel_pyr_slices_(std::shared_ptr<MChannelPyrSlices>(nullptr)),
                       file_format_(FileFormat::UNKNOWN), resolution_level_(-1),
                       xdim_(0), ydim_(0), zdim_(0), chunk_xdim_(0), chunk_ydim_(0), chunk_zdim_(0),
                       channel_pyr_dir_(std::string{}), slice_image_names_(std::unordered_map<std::string, std::vector<std::string>>{}),
                       voxel_type_(VoxelType::UNKNOWN), volume_complete_(VolumeComplete::UNKNOWN)  {}

    /// use channel_pyr_slices to construct new MChannelPyrSlices instance. do not share ownership of channel_pyr_slices
    MChannelPyrInfo(const mcp3d::MChannelPyrSlices &channel_pyr_slices):
            channel_pyr_slices_(std::make_shared<MChannelPyrSlices>(channel_pyr_slices)),
            file_format_(channel_pyr_slices_->file_format()), resolution_level_(-1),
            xdim_(0), ydim_(0), zdim_(0), chunk_xdim_(0), chunk_ydim_(0), chunk_zdim_(0),
            channel_pyr_dir_(channel_pyr_slices_->channel_pyr_dir()),
            slice_image_names_(std::unordered_map<std::string, std::vector<std::string>>{}),
            voxel_type_(VoxelType::UNKNOWN), volume_complete_(VolumeComplete::UNKNOWN)
    { ReadSliceImageNames(); }

    /// delegate to MChannelPyrInfo::MChannelPyrInfo(const MChannelPyrSlices& channel_pyr_slices)
    explicit MChannelPyrInfo(const std::string& channel_pyr_dir): MChannelPyrInfo(MChannelPyrSlices{channel_pyr_dir}) {}

    /// performs a deep copy of the MChannelPyrSlices object managed by other.channel_pyr_slices_
    /// do not share MChannelPyrSlices ownership with other
    MChannelPyrInfo(const MChannelPyrInfo& other);

    /// move constructor. channel_pyr_slices_ share ownership with objects previously
    /// sharing ownership with other.channel_pyr_slices_
    MChannelPyrInfo(MChannelPyrInfo&& other) noexcept;

    /// deep comparison of objects managed by channel_pyr_slices_
    /// does not compare the description_ string
    bool operator==(const MChannelPyrInfo& other) const;

    bool operator!=(const MChannelPyrInfo& other) const
    { return !(*this == other); }

    /// call channel_pry_slices.reset()
    void Clear();

    /// return true if channel_pyr_slices_ manages no object, or if channel_pyr_slices_ is empty
    bool empty() const;

    /// its possible for instance to have undefined dimensions. e.g. after instance is constructed with const string& and before
    /// MImageFormats::ReadChannelPyrInfo
    bool DimensionsDefined() const
    { return !empty() && zdim_ > 0 && ydim_ > 0 && xdim_ > 0; }

    /// return if not DimensionsDefined()
    /// check each slice. if any slice has less than expected number of images, volume is incomplete
    void UpdateVolumeCompleteness();

    /// return if empty()
    /// for each slice, read a vector of image names consistent with file_format_
    /// if file_format_ is unknown, create empty vector for the slice
    /// call UpdateVolumeCompleteness
    void ReadSliceImageNames();

    bool ZyxInVolume(int z, int y, int x) const
    { return z >= 0 && z < zdim() && y >= 0 && y < ydim() && x >= 0 && x < xdim(); }

    /// return number of entries present in slice_image_names_[slice_name],
    /// where slice_name is determined by MChannelPyrSlices::SliceName
    /// this is the number of image chunks actually written to the drive, not
    /// the number of image chunks expected in a fully written pyr level image
    int NumberOfImageChunksInSlice(int zslice_id, int yslice_id, int xslice_id) const;

    /// return number of image files expected in the slice if the pyr level image is fully written
    /// calculated from image chunk dimension, slice dimension and image dimension
    /// the edge slices (with highest slice id) can have fewer expected chunks than non edge slices
    int ExpectedNumberOfImageChunksInSlice(int zslice_id, int yslice_id, int xslice_id) const;

    int TotalSliceNumber() const
    { return slice_number(ChannelAxis::X) * slice_number(ChannelAxis::Y) * slice_number(ChannelAxis::Z); }

    /// obtain total number of file chunks by summing number of files under
    /// slice_image_names_. when instance is completely configured, this value
    /// should be equal to n_chunks_
    int TotalNumberOfImageChunks() const;

    long VolumeVoxels() const
    { return (long)xdim() * (long)ydim() * (long)zdim(); }

    long ChunkVoxels() const
    { return (long)chunk_xdim() * (long)chunk_ydim() * (long)chunk_zdim(); }

    long VolumeBytes() const
    { return VolumeVoxels() * (long) VoxelBytes(); }

    long ChunkBytes() const
    { return ChunkVoxels() * (long) VoxelBytes(); }

    int VoxelBytes() const
    { return empty() ? 0 : mcp3d::VoxelSize(voxel_type_); }

    int slice_number(ChannelAxis axis) const
    { return empty() ? 0 : channel_pyr_slices_->slice_number(axis); }

    /// if empty(), return 0
    /// if MChannelPyrSlices has dimension MAX_AXIS_SLICE_DIM for the axis, return dim(axis)
    int slice_dim(ChannelAxis axis) const;

    FileFormat file_format() const
    { return empty() ? mcp3d::FileFormat::UNKNOWN : file_format_; }

    int resolution_level() const
    { return empty() ? -1 : resolution_level_; }

    int zdim() const
    { return empty() ? 0: zdim_; }

    int ydim() const
    { return empty() ? 0 : ydim_; }

    int xdim() const
    { return empty() ? 0 : xdim_; }

    int dim(ChannelAxis axis) const;

    std::vector<int> xyz_dims() const
    { return std::vector<int>({zdim(), ydim(), xdim()}); }

    std::vector<int> xyz_chunk_dims() const
    { return std::vector<int>({chunk_zdim_, chunk_ydim_, chunk_xdim_}); }

    VoxelType voxel_type() const  { return voxel_type_; }

    int chunk_zdim() const
    { return empty() ? 0 : chunk_zdim_; }

    int chunk_ydim() const
    { return empty() ? 0 : chunk_ydim_; }

    int chunk_xdim() const
    { return empty() ? 0 : chunk_xdim_; }

    int chunk_dim(ChannelAxis axis) const;

    std::vector<int> chunk_xyz_dims() const
    { return std::vector<int>({chunk_zdim(), chunk_ydim(), chunk_xdim()}); }

    /// if zdim is not divisible by chunk_zdim, an extra chunk is needed along z axis to cover the residual zdim
    int n_zchunks() const
    { return DimensionsDefined() ? zdim_ / chunk_zdim_ + (int)(zdim_ % chunk_zdim_ > 0) : 0; }

    int n_ychunks() const
    { return DimensionsDefined() ? ydim_ / chunk_ydim_ + (int)(ydim_ % chunk_ydim_ > 0) : 0; }

    int n_xchunks() const
    { return DimensionsDefined() ? xdim_ / chunk_xdim_ + (int)(xdim_ % chunk_xdim_ > 0) : 0; }

    int n_chunks(ChannelAxis axis) const;

    int n_total_chunks()  const
    { return n_zchunks() * n_ychunks() * n_xchunks(); }

    const MChannelPyrSlices& channel_pyr_slices() const { return *(channel_pyr_slices_.get()); }

    std::string channel_pyr_dir() const  { return channel_pyr_dir_; };

    const std::unordered_map<std::string, std::vector<std::string>>& slice_image_names() const
    { return slice_image_names_; }

    VolumeComplete volume_complete() const  { return volume_complete_; }

    friend class MChannelInfo;
    friend class MImageFormats;
    friend class MTiffFormat;
    friend class MOmeTiffFormat;
    friend class MImarisFormat;
    friend class MHdf5Format;

private:
    /// copy construct channel_pyr_slices_ from channel_pyr_slices
    /// call ParseJson and AssertValid
    explicit MChannelPyrInfo(const nlohmann::json& info_json, const std::shared_ptr<MChannelPyrSlices>& channel_pyr_slices = std::shared_ptr<MChannelPyrSlices>(nullptr));

    const std::shared_ptr<MChannelPyrSlices>& channel_pyr_slices_ptr() const
    { return channel_pyr_slices_; }

    /// the coordinates (z, y, x) is local to the pyr_level
    /// assert DimensionsDefined()
    /// should only be called with volume_complete_ == VolumeComplete::COMPLETE (handling of incomplete volumes is through MImageFormats)
    /// return empty string for (z, y, x) outside of pyr_level volume range. (its possible for an in volume global coordinate to map to outside of
    /// pyr_level volume. use this function through MImageFormats for global volume range check)
    /// assuming the image names when sorted preserve zyx ordering of underlying data
    std::string ImagePath(int z, int y, int x) const;

    /// if a single slice cover the entire axis, return 1
    /// otherwise, if the slice_id is not the largest possible along the axis,
    /// return slice_dim / chunk dim. if slice_id is the largest possible, return
    /// the number of chunks required to cover the residual dim within the last slice
    int ExpectedNumberOfImageChunksInSliceAlongAxis(ChannelAxis axis, int slice_id) const;

    static FileFormat ParseFileFormat(const nlohmann::json& info_json);

    static int ParseResolutionLevel(const nlohmann::json& info_json);

    static int ParseDim(const nlohmann::json& info_json, ChannelAxis axis);

    static int ParseChunkDim(const nlohmann::json& info_json, ChannelAxis axis);

    /// make /a/b/c/ conisdered identical to /a/b/c, which is likely a user's intent
    static std::string ParseChannelPyrDir(const nlohmann::json& info_json);

    static VoxelType ParseVoxelType(const nlohmann::json& info_json);

    /// modify slice_image_names_ in place since it's potentially large
    void ParseSliceImageNames(const nlohmann::json& info_json);

    /// if channel_pyr_slices_ manages no object, construct from channel_pyr_dir_ parsed from json
    /// does not perform validity checks on parsed information. user should call AssertValid following call
    /// to ParseJson instead
    void ParseJson(const nlohmann::json& info_json);

    /// if empty() return
    /// (1) voxel_type_ is not unknown
    /// (2) channel_pyr_dir_ is equal to channel_pyr_slices_->channel_pyr_dir() and is valid directory
    /// (3) file_format_ is equal to channel_pyr_slices_->file_format()
    /// (4) resolution_level_ and channel_pyr_dir_ consistent (unless IMARIS format)
    /// if !DimensionDefined(), return
    /// (5) slice_dim * slice_number should be greater or equal to chunk_dim * n_chunks
    /// (6) dim, chunk dim and slice dim should all be positive.
    /// (7) dim should be greater than or equal to chunk dim and slice dim
    /// (8) slice dim should be greather than or equal to chunk dim, and divisible by chunk dim
    /// (9) chunk_dim * n_chunks can be greater than dim
    ///     e.g. xdim = 1100, chunk_xdim = 125, then 9 chunks along x axis are required
    ///     to cover xdim. the edge at the last xchunks are not actual data.
    ///     chunk_dim * (n_chunks - 1) should be strictly less than dim
    /// (10) slice_image_names_ key set is equal to channel_pyr_slices_->SliceNames() set
    /// (11) volume_complete_ is not unknown
    /// (10) expected number of image chunks in slice >=/== number of image chunks for the slice in slice_image_names_,
    ///      depending value of volume_complete_
    void AssertValid() const;

    nlohmann::json Jsonize() const;

    std::shared_ptr<MChannelPyrSlices> channel_pyr_slices_;
    FileFormat file_format_;
    int resolution_level_;
    int xdim_, ydim_, zdim_;
    int chunk_xdim_, chunk_ydim_, chunk_zdim_;
    std::string channel_pyr_dir_;
    /// slice path relative to channel pyramid level dir: image names under slice
    /// e.g. "z0000_000256/y0001_001024/x0002_001024":
    /// {"img_z0_y1024_x_2048", "img_z0_y1024_x_4096",
    ///  "img_z0_y1536_x_2048", "img_z0_y1536_x_4096"}
    /// if the channel pyramid directory's folder slice structure is not flat,
    /// each slice nust have the same number of files. padding files should be
    /// created when needed
    std::unordered_map<std::string, std::vector<std::string>> slice_image_names_;
    VoxelType voxel_type_;
    /// set to VolumeComplete::COMPLETE if all slices have expected number of image names
    /// if some slices are missing some image names, set to VolumeComplete::IMCOMPLETE
    VolumeComplete volume_complete_;
};

// xyz volume with pyramid hierarchy layouts
class MChannelInfo
{
public:
    MChannelInfo() = default;

    /// construct channel_layout_ from channel_root_dir. the instance will be considered the owner of channel_layout_
    explicit MChannelInfo(const std::string& channel_root_dir):
            channel_layout_(std::make_shared<MChannelLayout>(channel_root_dir)), channel_pyr_infos_(std::unordered_map<int, MChannelPyrInfo> {}) {}

    MChannelInfo(const MChannelInfo& other);

    MChannelInfo(MChannelInfo&& other) noexcept: channel_layout_(other.channel_layout_), channel_pyr_infos_(other.channel_pyr_infos_) {}

    /// return if pyr_info is empty
    /// call channel_layout_.ReadChannelLayout. assert pyr_info.resolution_level_
    /// exists for channel_layout_, and MChannelPyrSlices under pyr_info and channel_layout_ are equal
    /// if MChannelPyrInfo of same resolution level already exists, overwrite
    /// copy assign channel_layout_->channel_pyr_slices_ptr(pyr_level) to
    /// channel_pyr_infos_.at(pyr_level).channel_pyr_slices_
    /// pyr_info will be moved from
    void AddPyrInfo(MChannelPyrInfo &&pyr_info);

    // remove channel_pyr_infos_[pyr_level]. do nothing if pyr_level not in channel_pyr_infos_
    void RemovePyrInfo(int pyr_level)
    { channel_pyr_infos_.erase(pyr_level); }

    /// channel_layout_::ReadChannelLayout(). if a pyr_level no longer exists, remove
    /// it from channel_pyr_infos_
    void RefreshChannelLayout();

    bool operator==(const MChannelInfo& other) const
    { return *(channel_layout_.get()) == *(other.channel_layout_.get()) &&
             channel_pyr_infos_ == other.channel_pyr_infos_; }

    bool operator!=(const MChannelInfo& other) const
    { return !(*this == other); }

    bool HasPyrLevel(int pyr_level) const
    { return channel_layout_->HasPyrLevel(pyr_level); }

    bool HasPyrInfo(int pyr_level) const
    { return channel_pyr_infos_.find(pyr_level) != channel_pyr_infos_.end(); }

    bool LevelPathsExist(int pyr_level) const;

    /// high level sequence of writing operations should be such that images
    /// are written first, then image info. therefore no matter if MImage was
    /// constructed from storage, at image info saving all image paths must
    /// exist. image writing operation should update image info. image info
    /// writing is synchronized across threads and done only from rank 0 thread
    void Save() const;

    /// read saved __channel_info__.json
    /// for each pyr level, construct MChannelPyrInfo with the loaded json and channel_layout_->channel_pyr_slices(pyr_level)
    /// AssertValid() assuming the newly constructed MChannelPyrInfo from __channel_info__.json will be content of channel_pyr_infos_
    /// if no error, replace channel_pyr_infos_ with the newly constructed MChannelPyrInfo instances
    /// otherwise, leave instance in the same state as before Load()
    void Load();

    /// channel_pyr_infos_.clear()
    void ClearChannelPyrInfos()
    { channel_pyr_infos_.clear(); }

    int TotalSliceNumber(int pyr_level = 0) const
    { return HasPyrInfo(pyr_level) ? channel_pyr_info(pyr_level).TotalSliceNumber() : 0; }

    int TotalNumberOfImageChunks(int pyr_level = 0) const
    { return HasPyrInfo(pyr_level) ? channel_pyr_info(pyr_level).TotalNumberOfImageChunks() : 0; }

    long VolumeVoxels(int pyr_level = 0) const
    { return HasPyrInfo(pyr_level) ? channel_pyr_info(pyr_level).VolumeVoxels() : 0; }

    long ChunkVoxels(int pyr_level = 0) const
    { return HasPyrInfo(pyr_level) ? channel_pyr_info(pyr_level).ChunkVoxels() : 0; }

    long VolumeBytes(int pyr_level = 0) const
    { return HasPyrInfo(pyr_level) ? channel_pyr_info(pyr_level).VolumeBytes() : 0; }

    long ChunkBytes(int pyr_level = 0) const
    { return HasPyrInfo(pyr_level) ? channel_pyr_info(pyr_level).ChunkBytes() : 0; }

    int VoxelBytes(int pyr_level = 0) const
    { return HasPyrInfo(pyr_level) ? channel_pyr_info(pyr_level).VoxelBytes() : 0; }

    /// if output_pyr_level == 0, return {z, y, x}
    /// else if pyr_ratios()[output_pyr_level] does not exist, return empty vector
    /// does not check if z, y, x is in range or if output zyx is in range (even if zyx is negative)
    std::vector<int> ZyxGlobalToLocal(int output_pyr_level, int z, int y, int x) const
    { return ZyxTransform(output_pyr_level, z, y, x, true); }

    std::vector<int> ZyxLocalToGlobal(int input_pyr_level, int z, int y, int x) const
    { return ZyxTransform(input_pyr_level, z, y, x, false); }

    /// return channel_layout_->n_pyr_levels()
    int n_pyr_levels() const
    { return channel_layout_->n_pyr_levels(); }

    /// number of entries in channel_pyr_infos_
    int n_pyr_infos() const
    { return (int)channel_pyr_infos_.size(); }

    /// return sorted vector of pyr_levels in channel_pyr_infos_
    std::vector<int> pyr_info_levels() const;

    const MChannelLayout& channel_layout() const
    { return *(channel_layout_.get()); }

    const std::unordered_map<int, MChannelPyrInfo>& channel_pyr_infos() const
    { return channel_pyr_infos_; }

    const MChannelPyrInfo& channel_pyr_info(int pyr_level = 0) const;

    const MChannelPyrSlices& channel_pyr_slices(int pyr_level = 0) const;

    int z_scale_start_level() const;

    int slice_dim(int pyr_level, ChannelAxis axis) const
    { return HasPyrInfo(pyr_level) ? channel_pyr_info(pyr_level).slice_dim(axis) : 0; }

    int slice_number(int pyr_level, ChannelAxis axis) const
    { return HasPyrInfo(pyr_level) ? channel_pyr_info(pyr_level).slice_number(axis) : 0; }

    std::vector<int> xyz_dims(int pyr_level = 0) const
    { return HasPyrInfo(pyr_level) ? channel_pyr_info(pyr_level).xyz_dims() : std::vector<int>({0, 0, 0}); }

    std::vector<int> xyz_chunk_dims(int pyr_level = 0) const
    { return HasPyrInfo(pyr_level) ?  channel_pyr_info(pyr_level).xyz_chunk_dims() : std::vector<int>({0, 0, 0}); }

    int zdim(int pyr_level = 0) const
    { return xyz_dims(pyr_level)[0]; }

    int ydim(int pyr_level = 0) const
    { return xyz_dims(pyr_level)[1]; }

    int xdim(int pyr_level = 0) const
    { return xyz_dims(pyr_level)[2]; }

    int zchunk_dim(int pyr_level = 0) const
    { return xyz_chunk_dims(pyr_level)[0]; }

    int ychunk_dim(int pyr_level = 0) const
    { return xyz_chunk_dims(pyr_level)[1]; }

    int xchunk_dim(int pyr_level = 0) const
    { return xyz_chunk_dims(pyr_level)[2]; }

    FileFormat file_format(int pyr_level = 0) const
    { return HasPyrInfo(pyr_level) ? channel_pyr_info(pyr_level).file_format() : FileFormat::UNKNOWN; }

    VoxelType voxel_type(int pyr_level = 0) const
    { return HasPyrInfo(pyr_level) ? channel_pyr_info(pyr_level).voxel_type() : VoxelType::UNKNOWN; }

    std::string channel_root_dir() const
    { return channel_layout_->channel_root_dir(); }

    std::string channel_pyr_dir(int pyr_level = 0) const
    { return HasPyrInfo(pyr_level) ? channel_layout_->channel_pyr_dir(pyr_level) : std::string{}; }

    std::unordered_map<int, int> pyr_ratios(mcp3d::ChannelAxis axis) const;

    friend class MImageInfo;
    friend class MImageIO;

private:
    /// should only be called by MImageInfo instance. establish shared ownership of MChannelLayout instance between
    /// MImageInfo.volume_layout_ and this->channel_layout_
    explicit MChannelInfo(const std::shared_ptr<MChannelLayout>& channel_layout):
            channel_layout_(channel_layout), channel_pyr_infos_(std::unordered_map<int, MChannelPyrInfo> {}) {}

    /// does not assert returned path exists
    std::string channel_info_path() const
    { return mcp3d::JoinPath(channel_layout_->channel_root_dir_, mcp3d::MChannelInfo::channel_info_json_name()); }

    // __channel_info__.json
    static std::string channel_info_json_name()
    { return "__channel_info__.json"; }

    /// if empty(), return
    /// if n_pyr_level() == 0, return
    /// (1) call channel_layout_->AssertValid
    /// (2) assert all pyr levels in pyr_info_levels() exist in channel_layout_
    ///     (its possible to have pyr_level in channel_layout_ that does not exist in channel_pyr_infos_)
    /// (3) pyramid ratios must be geometric series of 2^i, i = {0, 1, 2...}
    ///     for the z axis, allow the scaling to occur deeper into the pyrmaid
    ///     levels. additionally, allow the ratio to be equal to the preceding
    ///     scaled level's zdim, as the z axis often has smaller range compared
    ///     to xy axis and would degenerate towards zero if the geometric series
    ///     is strictly reinforced
    ///     e.g. zdims of 512, 512, 256, 128, 64, 32, 32, 32 is valid. the same
    ///     values are not valid for xy axis.
    ///     pyramid ratio at each level should be no greater than the dimension along any axis
    ///     pyramid ratios will only be validated if level 0 MChannelPyrInfo exists
    ///     for the instance (pyr_ratios can not be calculated otherwise)
    /// (4) call AssertValid on each MChannelPyrInfo instance and assert non empty
    /// (5) assert that channel_layout_ and channel_pyr_infos_ shared_ptr stores same pointers to MChannelPyrSlices instance
    void AssertValid() const;

    std::vector<int> ZyxTransform(int pyr_level, int z, int y, int x, bool global_to_local) const;

    /// if not HasPyrInfo(0), do nothing.
    /// if at a given pyr_level, the pyr_ratio along an axis is greater than the dimension of the axis,
    /// remove channel_pyr_infos_[pyr_level]
    void RemovePyrInfoWithUndersizedVolumes();

    const std::shared_ptr<MChannelLayout>& channel_layout_ptr() const
    { return channel_layout_; }

    std::shared_ptr<MChannelLayout> channel_layout_;
    // pyr_level: MChannelPyrInfo
    std::unordered_map<int, MChannelPyrInfo> channel_pyr_infos_;
};
}

#endif //MCP3D_MCP3D_CHANNEL_INFO_HPP
