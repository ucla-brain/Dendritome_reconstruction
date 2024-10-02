//
// Created by muyezhu on 2/26/18.
//
#include <algorithm>
#include "image_layout/mcp3d_channel_pyr_slices.hpp"
#include "mcp3d_image_utils.hpp"
#include "mcp3d_image_view.hpp"

using namespace std;

mcp3d::MImageBlock::MImageBlock(const vector<int>& offsets, const vector<int>& extents, const vector<int>& strides)
{
    MCP3D_ASSERT(offsets.size() <= (size_t)3 && extents.size() <= (size_t)3 && strides.size() <= (size_t)3)
    if (! offsets.empty() && ! AllNonNegative(offsets))
        MCP3D_OUT_OF_RANGE("negative offsets encountered")
    if (! extents.empty() && ! AllNonNegative(extents))
        MCP3D_OUT_OF_RANGE("negative extents encountered")
    if (! strides.empty() && !AllPositive(strides))
        MCP3D_OUT_OF_RANGE("non positive strides encountered")

    offsets_ = make_unique<int[]>(3);
    for (size_t i = 0; i < 3; ++i)
    {
        size_t offset_pos = i + offsets.size() - 3;
        if (offset_pos >= 0 and offset_pos < offsets.size())
            offsets_[i] = offsets[offset_pos];
        else
            offsets_[i] = 0;
    }
    extents_ = make_unique<int[]>(3);
    for (size_t i = 0; i < 3; ++i)
    {
        size_t extent_pos = i + extents.size() - 3;
        if (extent_pos >= 0 and extent_pos < extents.size())
            extents_[i] = extents[extent_pos];
        else
            extents_[i] = 0;
    }
    strides_ = make_unique<int[]>(3);
    for (size_t i = 0; i < 3; ++i)
    {
        size_t stride_pos = i + strides.size() - 3;
        if (stride_pos >= 0 and stride_pos < strides.size())
            strides_[i] = strides[stride_pos];
        else
            strides_[i] = 1;
    }
}

mcp3d::MImageBlock::MImageBlock(const mcp3d::MImageBlock& other)
{
    offsets_ = make_unique<int[]>(3);
    extents_ = make_unique<int[]>(3);
    strides_ = make_unique<int[]>(3);
    for (int i = 0; i < 3; ++i)
    {
        offsets_[i] = other.offsets_[i];
        extents_[i] = other.extents_[i];
        strides_[i] = other.strides_[i];
    }
}

mcp3d::MImageBlock& mcp3d::MImageBlock::operator=(const mcp3d::MImageBlock& other)
{
    for (int i = 0; i < 3; ++i)
    {
        offsets_[i] = other.offsets_[i];
        extents_[i] = other.extents_[i];
        strides_[i] = other.strides_[i];
    }
    return *this;
}

bool mcp3d::MImageBlock::operator== (const mcp3d::MImageBlock& other) const
{
    for (int i = 0; i < 3; ++i)
        if (offsets_[i] != other.offsets_[i])
            return false;
    for (int i = 0; i < 3; ++i)
        if (extents_[i] != other.extents_[i])
            return false;
    for (int i = 0; i < 3; ++i)
        if (strides_[i] != other.strides_[i])
            return false;
    return true;
}

void mcp3d::MImageBlock::Clear()
{
    for (int i = 0; i < 3; ++i)
    {
        offsets_[i] = 0;
        extents_[i] = 0;
        strides_[i] = 1;
    }
}

void mcp3d::MImageBlock::PrintBlock() const
{
    cout << "offsets(zyx): " << mcp3d::JoinVector(offsets(), ", ", true) << ", ";
    cout << "extents(zyx): " << mcp3d::JoinVector(extents(), ", ", true) << ", ";
    cout << "strides(zyx): " << mcp3d::JoinVector(strides(), ", ", true) << endl;
}

mcp3d::MImageView::MImageView(const MImageView &other): image_info_(other.image_info_)
{
    pyr_level_offsets_ = other.pyr_level_offsets_;
    pyr_level_extents_ = other.pyr_level_extents_;
    pyr_level_strides_ = other.pyr_level_strides_;
    global_image_block_ = other.global_image_block_;
    view_channels_ = other.view_channels_;
    view_times_ = other.view_times_;
    pyr_level_ = other.pyr_level_;
    interpret_block_as_local_ = other.interpret_block_as_local_;
    voxel_type_ = other.voxel_type_;
}

bool mcp3d::MImageView::operator==(const mcp3d::MImageView &other) const
{
    return image_info_ == other.image_info_ &&
           pyr_level_offsets_ == other.pyr_level_offsets_ &&
           pyr_level_extents_ == other.pyr_level_extents_ &&
           pyr_level_strides_ == other.pyr_level_strides_ &&
           view_channels_ == other.view_channels_ &&
           view_times_ == other.view_times_ &&
           pyr_level_ == other.pyr_level() &&
           voxel_type_ == other.voxel_type() &&
           interpret_block_as_local_ == other.interpret_block_as_local_;
}

void mcp3d::MImageView::SelectView(const MImageBlock &image_block, const vector<int>& channels, int pyr_level, bool interpret_block_as_local)
{
    // clear view
    Clear();
    try
    {
        MCP3D_ASSERT(!image_info_.empty())
        // update member data fields
        pyr_level_ = pyr_level;
        view_channels_ = channels;
        MCP3D_ASSERT(!view_channels_.empty())
        for (int view_channel: view_channels_)
        {
            if (!(view_channel >= 0 && view_channel < image_info_.volume_layout().n_channels()))
                MCP3D_RUNTIME_ERROR("channel number " + to_string(view_channel) + " is out of range")
            if (!image_info_.HasChannelPyrInfo(view_channel, 0))
                MCP3D_RUNTIME_ERROR("no MChannelPyrInfo exists at pyr_level 0 for channel " + to_string(view_channel))
            if (!image_info_.HasChannelPyrInfo(view_channel, pyr_level_))
                MCP3D_RUNTIME_ERROR("no MChannelPyrInfo exists at pyr_level " + to_string(pyr_level_) + " for channel " + to_string(view_channel))
        }
        view_times_ = vector<int>(1, 0);
        interpret_block_as_local_ = interpret_block_as_local;

        // since MImageInfo::AssertValid would have ensured MChannelInfo instances to be compatible, can use
        // first view channel to retrieve MChannelPyrInfo properties
        int channel = view_channels_[0];
        // short alias
        const mcp3d::MChannelInfo& channel_info = image_info_.channel_info(channel, 0);
        // if voxel type unknown, use storage voxel type of first view channel
        voxel_type_ = channel_info.voxel_type(pyr_level);

        pyr_level_offsets_ = vector<int>(3, 0);
        pyr_level_extents_ = vector<int>(3, 0);
        pyr_level_strides_ = vector<int>(3, 1);

        if (interpret_block_as_local)  // image block is at pyr_level
        {
            for (int i = 0; i < 3; ++i)
            {
                pyr_level_offsets_[i] = image_block.offsets_[i];
                if (image_block.extents_[i] == 0)
                    pyr_level_extents_[i] = channel_info.xyz_dims(pyr_level_)[i] - pyr_level_offsets_[i];
                else
                    pyr_level_extents_[i] = image_block.extents_[i];
                pyr_level_strides_[i] = image_block.strides_[i];
            }
            for (int i = 0; i < 3; ++i)
            {
                int upscale_factor = (i == 0) ? channel_info.pyr_ratios(mcp3d::ChannelAxis::Z)[pyr_level_] :
                                     (i == 1) ? channel_info.pyr_ratios(mcp3d::ChannelAxis::Y)[pyr_level_] : channel_info.pyr_ratios(mcp3d::ChannelAxis::X)[pyr_level_];
                global_image_block_.offsets_[i] = pyr_level_offsets_[i] * upscale_factor;
                global_image_block_.extents_[i] = pyr_level_extents_[i] * upscale_factor;
                // not scaling up strides unless strides is greater than one.
                // there's multiple possible level 0 strides that can result in
                // current level stride equal to 1
                global_image_block_.strides_[i] = pyr_level_strides_[i] == 1 ? 1 : pyr_level_strides_[i] * upscale_factor;
            }
        }
        else // image block is at level 0
        {
            for (int i = 0; i < 3; ++i)
            {
                global_image_block_.offsets_[i] = image_block.offsets_[i];
                if (image_block.extents_[i] == 0)
                    global_image_block_.extents_[i] = channel_info.xyz_dims(0)[i] - global_image_block_.offsets_[i];
                else
                    global_image_block_.extents_[i] = image_block.extents_[i];
                global_image_block_.strides_[i] = image_block.strides_[i];
            }

            for (int i = 0; i < 3; ++i)
            {
                int downscale_factor = i == 0 ? channel_info.pyr_ratios(mcp3d::ChannelAxis::Z)[pyr_level_] :
                                       (i == 1) ? channel_info.pyr_ratios(mcp3d::ChannelAxis::Y)[pyr_level_] : channel_info.pyr_ratios(mcp3d::ChannelAxis::X)[pyr_level_];
                pyr_level_offsets_[i] = global_image_block_.offsets_[i] / downscale_factor;
                pyr_level_extents_[i] = max(1, global_image_block_.extents_[i] / downscale_factor);
                pyr_level_strides_[i] = max(1, global_image_block_.strides_[i] / downscale_factor);
            }
        }

        // range validation for global block for xyz
        for (int i = 0; i < 3; ++i)
        {
            int level0_dim = channel_info.xyz_dims(0)[i];
            if (global_image_block_.offsets_ptr()[i] >= level0_dim)
                MCP3D_OUT_OF_RANGE("global image block: index " + to_string(global_image_block_.offsets_ptr()[i]) +
                                   " is out of bounds for axis " + to_string(i) + " with size " + to_string(level0_dim))
        }
    }
    catch (...)
    {
        Clear();
        MCP3D_RETHROW(current_exception())
    }
}

bool mcp3d::MImageView::OutOfPyrImageBoundary() const
{
    if (empty())
        return false;
    int view_channel = view_channels_[0];
    return !mcp3d::AllGreater(image_info_.channel_info(view_channel).xyz_dims(pyr_level_), pyr_level_offsets_);
}

bool mcp3d::MImageView::PartiallyOutOfPyrImageBoundary() const
{
    if (empty())
        return false;
    // if view entirely out of boundary, return true
    if (OutOfPyrImageBoundary())
        return true;
    int view_channel = view_channels_[0];
    vector<int> n_in_boundary_voxels = mcp3d::StridedExtents(mcp3d::SubtractSeq<int>(image_info_.channel_info(view_channel, 0).xyz_dims(pyr_level_), pyr_level_offsets_),
                                                             pyr_level_strides_);
    vector<int> in_boundary_full_stride_extents = (n_in_boundary_voxels - 1) * pyr_level_strides_ + 1;
    vector<int> remaining_extents = pyr_level_extents_ - in_boundary_full_stride_extents;
    return !mcp3d::AllGreater(pyr_level_strides_, remaining_extents);
}


void mcp3d::MImageView::Clear()
{
    global_image_block_.Clear();
    pyr_level_offsets_.clear();
    pyr_level_extents_.clear();
    pyr_level_strides_.clear();
    view_channels_.clear();
    view_times_.clear();
    pyr_level_ = -1;
    voxel_type_ = mcp3d::VoxelType::UNKNOWN;
}

void mcp3d::MImageView::PrintView() const
{
    if (global_image_block_.extents_[0] == 0)
    {
        cout << "global image view uninitialized" << endl;
        return;
    }
    cout << "global image volume:" << endl;
    global_image_block_.PrintBlock();
    cout << "selected view: pyramid level = " << pyr_level_ << ", "
         << "time points = " << mcp3d::JoinVector<int>(view_times(), ", ") << ", "
         << "channels = " << mcp3d::JoinVector<int>(view_channels(), ", ") << endl;
    cout << "offsets(zyx) = ";
    cout << mcp3d::JoinVector(pyr_level_offsets_, ", ", true);
    cout << ", ";
    cout << "extents(zyx) = ";
    cout << mcp3d::JoinVector(pyr_level_extents_, ", ", true);
    cout << ", ";
    cout << "strides(zyx) = ";
    cout << mcp3d::JoinVector(pyr_level_strides_, ", ", true);
    cout << endl;
    cout << "equal to " << ViewMemorySize("GB") << " GB of memory" << endl;
}

bool mcp3d::MImageView::is_unit_strided(const string &axes) const
{
    if (empty())
        return false;
    int strides = 1;
    string axes_ = mcp3d::StringLower(axes);
    for (const auto& axis: axes_)
    {
        if (axis == 'x')
            strides *= pyr_level_strides_[2];
        else if (axis == 'y')
            strides *= pyr_level_strides_[1];
        else if (axis == 'z')
            strides *= pyr_level_strides_[0];
    }
    return strides == 1;
}

vector<int> mcp3d::MImageView::xyz_dims() const
{ return empty() ? std::vector<int>({0, 0, 0}) : mcp3d::StridedExtents(pyr_level_extents_, pyr_level_strides_); }

vector<int> mcp3d::MImageView::xyz_dims_in_bound() const
{
    if (empty() || OutOfPyrImageBoundary())
        return vector<int>({0, 0, 0});
    vector<int> view_end = pyr_level_offsets_ + pyr_level_extents_;
    vector<int> view_valid_end = mcp3d::Minimum(view_end, image_info_.xyz_dims(pyr_level_));
    vector<int> view_in_bound_extents =
            mcp3d::StridedExtents(mcp3d::SubtractSeq<int>(view_valid_end, pyr_level_offsets_), pyr_level_strides_);
    return view_in_bound_extents;
}

vector<int> mcp3d::MImageView::dims() const
{
    if (empty())
        return {0, 0, 0, 0, 0};
    vector<int> xyz_dimensions(xyz_dims());
    return {1, n_channels(), xyz_dimensions[0], xyz_dimensions[1], xyz_dimensions[2]};
}

int mcp3d::MImageView::view_volume_index(int channel_number, bool throw_non_exist) const
{
    for (int i = 0; i < (int)view_channels_.size(); ++i)
        if (view_channels_[i] == channel_number)
            return i;
    if (throw_non_exist)
        MCP3D_RUNTIME_ERROR("channel_number " + to_string(channel_number) + " is not in view channels")
    return -1;
}

