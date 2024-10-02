//
// Created by muyezhu on 2/26/18.
//

#ifndef MCP3D_MCP3D_IMAGE_REGION_HPP
#define MCP3D_MCP3D_IMAGE_REGION_HPP

#include "common/mcp3d_common.hpp"
#include "mcp3d_image_info.hpp"

namespace mcp3d
{

class MImageBlock
{
public:
    /// offsets_, extents_ and strides_ contain 3 dimensions, in order of zyx
    /// if values missing from constructor parameter, offsets_ and extents_
    /// prepend 0, strides_ prepend 1.
    /// channel_ is same length as input when input
    /// is not empty. otherwise, if any dimension of extents equal to 0,
    /// n_channel_ is set to zero, else n_channel_ = 1 and channels[0] = 0
    /// times_ is a single 0 at the moment if n_channels_ > 0, else n_times_ = 0
    /// channels_ will always hold max_channel_number_ integers for correct code
    /// behavior if introducing more channels from copied instance or within
    /// view selection. similarly times_ will always hold max_time_numer_
    explicit MImageBlock(const std::vector<int>& offsets = std::vector<int>(),
                         const std::vector<int>& extents = std::vector<int>(),
                         const std::vector<int>& strides = std::vector<int>());

    MImageBlock(const MImageBlock& other);

    MImageBlock& operator=(const MImageBlock& other);

    bool operator== (const MImageBlock& other) const;

    bool empty() const  { return AnyZeros(extents()); }

    void Clear();

    std::vector<int> offsets() const  { return std::vector<int>({offsets_[0], offsets_[1], offsets_[2]}); }

    std::vector<int> extents() const  { return std::vector<int>({extents_[0], extents_[1], extents_[2]}); }

    std::vector<int> strides() const  { return std::vector<int>({strides_[0], strides_[1], strides_[2]}); }

    int* offsets_ptr()  { return offsets_.get(); }

    int* extents_ptr()  { return extents_.get(); }

    int* strides_ptr()  { return strides_.get(); }

    void PrintBlock() const;

    friend class MImageView;

private:
    std::unique_ptr<int[]> offsets_, extents_, strides_;
};

class MImageView
{
    /// manages the view into an image volume through the image pyramids
    /// offsets, extents and strides apply to zyx dimensions
    /// channels and times are vectors of selected channels and time points,
    /// does not need to be continuous in value
    /// view volumes across channels and time points have the same pyr level,
    /// offsets, extents and strides
public:
    explicit MImageView(): image_info_{MImageInfo{}}, global_image_block_{}, voxel_type_(VoxelType::UNKNOWN), pyr_level_(-1), interpret_block_as_local_(false),
                           pyr_level_offsets_(std::vector<int>{}), pyr_level_extents_(std::vector<int>{}), pyr_level_strides_(std::vector<int>{}),
                           view_channels_{}, view_times_{} {};

    explicit MImageView(const MImageInfo& image_info):
            image_info_{image_info}, global_image_block_{}, voxel_type_(VoxelType::UNKNOWN), pyr_level_(-1), interpret_block_as_local_(false),
            pyr_level_offsets_(std::vector<int>{}), pyr_level_extents_(std::vector<int>{}), pyr_level_strides_(std::vector<int>{}),
            view_channels_{}, view_times_{}
    { MCP3D_ASSERT(!image_info_.empty()) };

    MImageView(const MImageView& other);

    /// includes deep comparison of image_info_ and other.image_info_
    bool operator== (const MImageView& other) const;

    /// image_info_ can not be reassigned
    MImageView& operator=(const MImageView& other) = delete;

    MImageView& operator=(MImageView&& other) = delete;

    bool empty() const { return view_channels_.empty(); }

    /// const MImageBlock& image_block provides offsets, extents, strides in the global or local image
    /// coordinate system, depending on value of interpret_view_as_local. the offsets and extents define
    /// a rectangular volume within which voxels wll be retrieved. the strides parameter further
    /// determines if the voxels will be retrieved in strided manner.
    /// all channels and time points in view will have identical global offsets, extents, and strides.
    /// if the view is global, MImageView fills in default values for image_block if needed, and copy construct
    /// global_image_block_ data member from image_block. MImageView's own pyr_level_offsets_, pyr_level_extents_,
    /// pyr_level_strides_ are then calculated from the global coordinate system at requested pyramid level by scaling
    /// xyz according to pyramid ratio of the pyramid level.
    /// therefore retrieve image at given pyramid level with pyr_level_offsets_, pyr_level_extents_, pyr_level_strides_
    /// refers to the same level 0 image area identified by global_image_block_.
    /// if image_block has local pyramid level interpretation, global_image_block_ is produced
    /// from image_block in similar but reverse process
    /// in both global view and pyramid level view, default values for all offsets along all axes is 0
    /// default values for all extents along all axes is from offset till
    /// last element along the given axis. offsets + extent can exceed axis
    /// length. extents with zero value will be changed to default since empty
    /// view selection is not meaningful. the out of range portion of data will
    /// be padded with background pixels. but value of offset in global image block
    /// must be within axis range (at other view level the offsets are allowed to
    /// run out of bounds due to earlier level boundary voxel may not have correspondence in later levels)
    /// default values for all strides along all axes is 1
    /// channels can not be empty. image_info_ must have at least MChannelPyrInfo at levels pyr_level and 0 for the given channels
    /// calls Clear() immediately. if error is thrown later in function body, catch and call Clear(),
    /// then rethrow. therefore failed selection should result in empty() = true
    void SelectView(const MImageBlock &image_block, const std::vector<int>& channels, int pyr_level = 0, bool interpret_block_as_local = false);

    void SelectView(const MImageBlock &image_block, int channel, int pyr_level = 0, bool interpret_block_as_local = false)
    { SelectView(image_block, std::vector<int>({channel}), pyr_level, interpret_block_as_local); }

    void SelectView(const MImageView& other)
    { SelectView(other.global_image_block_, other.view_channels_, other.pyr_level_, other.interpret_block_as_local_); }

    /// if view is empty, return false
    /// if view is entirely out of pyramid image boundary return true
    bool OutOfPyrImageBoundary() const;

    /// if view is empty, return false
    /// if view is partially out of pyramid image boundary return true
    /// if ViewEntirelyOutOfPyrImageBoundary() is true, return true
    bool PartiallyOutOfPyrImageBoundary() const;

    /// clear view selection.
    void Clear();

    void PrintView() const;

    double ViewMemorySize(const std::string &unit) const
    { return empty() ? 0.0 : MemorySize(VolumeBytes(), unit) * (double)n_channels(); }

    const MImageBlock& global_image_block() const
    { return global_image_block_; }

    int pyr_level() const
    { return pyr_level_; }

    bool interpret_block_as_local() const
    { return interpret_block_as_local_; }

    std::vector<int> pyr_level_offsets() const
    { return pyr_level_offsets_; }

    /// the extents covered by the selected block at view level image,
    /// may not be the dimension of retrieved data due to striding, equal to
    /// xyz_dims() if strides are all 1
    std::vector<int> pyr_level_extents() const
    { return pyr_level_extents_; }

    /// the stride values used to read view level image
    std::vector<int> pyr_level_strides() const
    { return pyr_level_strides_; }

    std::vector<int> view_channels() const { return view_channels_; }

    // given a channel, return its index in view_channels_
    // for channel = 3 and view_channels_ = [0, 1, 3], return 2
    int view_volume_index(int channel_number, bool throw_non_exist = false) const;

    std::vector<int> view_times() const { return view_times_; }

    /// return true if unit stride along all queried axis. for empty() instance
    /// return false
    bool is_unit_strided(const std::string &axes = "xyz") const;

    int n_channels() const
    { return empty() ? 0 : (int)view_channels_.size(); }

    int n_times() const
    { return empty() ? 0 : (int)view_times_.size(); }

    /// number of xyz volumes in view
    int n_volumes() const
    { return n_channels() * n_times(); }

    /// if empty(), return {0, 0, 0}
    /// xyz dimensions of view data, accounting for striding
    /// define in cpp: mcp3d_image_utils.hpp has include from mcp3d_image_view.hpp
    std::vector<int> xyz_dims() const;

    /// xyz dimensions of view data that falls within the volume bounds at view level, accounting for striding
    /// if empty(), return {0, 0, 0}
    /// its possible to return all 0s vector from non empty instance: at volume boundaries there are voxels at level 0 that
    /// does not have a mapping at coarser levels. this is the dimensions of actual data retrieval from disk storage
    std::vector<int> xyz_dims_in_bound() const;

    /// tcxyz dimension of view data, accounting for striding
    std::vector<int> dims() const;

    int zdim() const
    { return xyz_dims()[0]; }

    int ydim() const
    { return xyz_dims()[1]; }

    int xdim() const
    { return xyz_dims()[2]; }

    /// number of voxels in an xyz volume
    long VolumeVoxels() const
    { return empty() ? 0 : (long)xdim() * (long)ydim() * (long)zdim(); }

    /// number of voxels in an xy plane
    long PlaneVoxels() const
    { return empty() ? 0 : (long)xdim() * (long)ydim(); }

    int VoxelBytes() const
    { return empty() ? 0 : VoxelSize(voxel_type_); }

    long VolumeBytes() const
    { return VolumeVoxels() * VoxelBytes(); }

    long PlaneBytes() const
    { return PlaneVoxels() * VoxelBytes(); }

    VoxelType voxel_type() const
    { return empty() ? VoxelType::UNKNOWN : voxel_type_; }

    friend class MImageBase;
    friend class MImage;
    friend class MImageInMemory;
    friend class MImageIO;

private:
    const MImageInfo& image_info_;
    MImageBlock global_image_block_;
    VoxelType voxel_type_;
    int pyr_level_;
    bool interpret_block_as_local_;
    // offsets, extents and strides vectors of view volumes
    std::vector<int> pyr_level_offsets_, pyr_level_extents_, pyr_level_strides_;
    // view_channels_: channels in view. for an image with 5 channels,
    // view_channels_ can be {0, 1, 4} etc.
    std::vector<int> view_channels_, view_times_;
};

}

#endif //MCP3D_MCP3D_IMAGE_REGION_HPP
