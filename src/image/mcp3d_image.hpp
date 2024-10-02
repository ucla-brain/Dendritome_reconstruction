//
// Created by muyezhu on 9/15/17.
// currently only plan to support tif
// if rgb channels encountered will convert to gray
//

#ifndef MCP3D_MCP3D_IMAGE_HPP
#define MCP3D_MCP3D_IMAGE_HPP

#include <iostream>
#include <memory>
#include <new>
#include <utility>
#include <vector>
#include <omp.h>
#include "common/mcp3d_common.hpp"
#include "image_layout/mcp3d_volume_layout.hpp"
#include "mcp3d_image_utils.hpp"
#include "mcp3d_image_info.hpp"
#include "mcp3d_image_view.hpp"
#include "mcp3d_image_base.hpp"


namespace mcp3d
{
class MImage: public MImageBase
{
public:
    MImage() = default;
    
    explicit MImage(const std::string &volume_path, const std::vector<std::string>& channel_dir_names = std::vector<std::string>{}):
            MImageBase{}, image_info_(volume_path, channel_dir_names), selected_view_{image_info_}, loaded_view_{image_info_} {}

    explicit MImage(MVolumeLayout&& volume_layout):
            MImageBase{}, image_info_(std::move(volume_layout)), selected_view_{image_info_}, loaded_view_{image_info_} {}

    explicit MImage(const MImageInfo& image_info): MImageBase{}, image_info_{image_info}, selected_view_{image_info_}, loaded_view_{image_info_} {}

    /// make a copy of the loaded data of other
    MImage(const MImage& other);

    bool empty() const override
    { return data_.empty(); }

    /// implementation at MImageIO will remove requested MChannelPyrInfo at pyr levels before loading/reading them
    /// if load_from_json is false, MChannelInfo will only be read from disk image files, previously saved
    /// channel info json files are ignored. this can be helpful when there exists multiple symbolic links
    /// pointing to the data directory. if pyr_level is negative, MChannelPyrInfo from all pyr_levels are
    /// read. otherwise, only MChannelInfo at pyr_level is read
    void ReadImageInfo(const std::vector<int> &channel_numbers, int pyr_level = -1, bool load_from_json = true);

    // convenient overload for reading single channel
    // note that if channel_number is greater than 0, the MVolumeLayout instance
    // must succeed to parse multiple channels
    void ReadImageInfo(int channel_number, int pyr_level = -1, bool load_from_json = false)
    {  ReadImageInfo(std::vector<int>({channel_number}), pyr_level, load_from_json); }

    /// if channel_numbers is an empty vector, attempt to select view for all channels present in volume_layout_
    /// if the MChannelPyrInfo at the specified channel numbers and pyr level has not been read, read it first
    /// if the pyr level is positive (> 0) and level 0 MChannelPyrInfo at the specified channel number has not
    /// been read, read level 0 MChannelPyrInfo as well
    /// select a view into the stored image data. the selected view is the
    /// staged view for which stored data is retrieved in next ReadData call
    /// the voxel type of the selected view is identical as stored voxel type.
    /// if different voxel type is desired, explicit casting should be performed
    /// the image view selected has global context
    /// the function gaurantees the offsets of the global image block is within
    /// level 0 image boundary. at other image pyramid level, the offsets are
    /// allowed to run out of boundary, in which case reading image will return
    /// all background voxels. this is due to an level 0 with for example width
    /// 11, where its level 0 image will have width 5. global view offset = 10
    /// along width is valid, but will be out of bounds (10 / 2 = 5) at level 1
    /// if allocate is true, call AllocateSelection. if MImageView instance sees
    /// empty channel_numbers vector, it attempts to select block from all channels
    /// under image_info_.volume_layout_
    void SelectView(const MImageBlock &view_block, const std::vector<int>& channel_numbers = std::vector<int> {},
                    int pyr_level = 0, bool interpret_block_as_local = false, bool allocate = false);

    /// convenient signature for single channel operations
    void SelectView(const MImageBlock &view_block, int channel_number, int pyr_level = 0,
                    bool interpret_block_as_local = false, bool allocate = false)
    { SelectView(view_block, std::vector<int>({channel_number}), pyr_level, interpret_block_as_local, allocate); }

    /// views partially outside of global image volume is valid. the out of volume
    /// voxels will be filled with valid. views have no overlap with global image volume are invalid
    /// if reading operation fails, the data vector is cleared since its now in
    /// unknown state. image loaded view is cleared. caught error is rethrown
    /// data is read in storage voxel type
    void ReadData(const std::string &mode = "verbose");

    /// clear data_ vector and loaded_view_
    void ClearLoaded() override;

    /// assert !empty() and volume index valid
    /// returns pointer stored in data_
    /// use channel number (not volume index) to retrieve data. channel number
    /// is converted to volume index internally. channel number is determined by
    /// volume_layout_ instance. volume_index is an index into the data_ vector
    template <typename VType = uint8_t>
    VType* Volume(int channel_number = 0, int time = 0);

    /// read only exposure to data pointer
    template <typename VType = uint8_t>
    const VType* ConstVolume(int channel_number = 0, int time = 0) const;

    template <typename VType = uint8_t>
    VType* Plane(int channel_number, int view_z, int time = 0);

    /// the coordinates here are local to voxels in loaded image view.
    /// first element in first array of the data_ vector is indexed (0, 0, 0, 0, 0)
    template <typename VType = uint8_t>
    void SetVoxel(int c, int z, int y, int x, VType val, int t = 0)
    { Volume<VType>(c, t)[VoxelAddress(c, z, y, x)] = val; }

    template <typename VType = uint8_t>
    VType& operator() (int c, int z, int y, int x, int t = 0)
    { return Volume<VType>(c, t)[VoxelAddress(c, z, y, x)]; }

    template <typename VType = uint8_t>
    const VType& At(int c, int z, int y, int x, int t = 0) const
    { return ConstVolume<VType>(c, t)[VoxelAddress(c, z, y, x)]; }

    /// should only be called by rank 0 process or thread.
    /// only the channels for which MImageChannelInfo has been read are saved
    void SaveImageInfo();

    /// write output with opencv
    void WriteViewXYPlane(const std::string &img_path, int c, int z, int t = 0);

    /// write loaded data. if no data loaded, load selected data first. if no
    /// selection is made either, do nothing
    /// output name should contain x, y, z start position of the view, with xyz
    /// values relative to the dimensions at the view's pyr_level. additionally
    /// stride, channel, time should also be in name
    void WriteViewVolume(const std::string &out_dir, const std::string& img_name_prefix, FileFormat volume_format = FileFormat::UNKNOWN);

    /// if this function called on local machine with remote host option, validate
    /// call image states, and qsub an MPI executable on remote. local
    /// execution will do multi threading unless instructed not to
    /// writing image pyramid successively, using start_parent_level as initial
    /// image volume to downsize from, while the new pyramid level to be written
    /// is less than end_level. if all is well the pyramid levels should be
    /// generated with [start_parent_level, end_parent_level) as parent image,
    /// leaving the image with pyramid levels
    /// [0, ..., start_parent_level, start_parent_level + 1, ..., end_parent_level + 1)
    /// if end_parent_level argument not given, all possible levels are written
    /// if image_info_.MaxMappablePyrLevel() is positive (can be calculated)
    /// rgb tiff image's pyramids will be gray
    void WriteImagePyramids(int channel_number, int start_parent_level = 0, int end_parent_level = 999,
                            bool multi_threading = false, bool save_image_info = true, FileFormat write_format = FileFormat::UNKNOWN);

    void WriteImagePyramid(int channel_number, int parent_level = 0, bool multi_threading = false,
                           bool save_image_info = true, FileFormat write_format = FileFormat::UNKNOWN);

    // acquire ownership of data in other_data and release current existing data
    // use with extreme caution: consistency between other_data and image_info_,
    // image view objects are not validated. mostly used for performance
    // considerations in certain operations
    template <typename VType = uint8_t>
    void AcquireData(std::vector<std::unique_ptr<VType[]>> &other_data);

    // volume_index = -1 will append other data to data_. if volume_index in
    // [0, n_volumes()), replace data at indicated position
    template <typename VType = uint8_t>
    void AcquireData(std::unique_ptr<VType[]> &other_data, int volume_index = -1);

    // release ownership of data. clear loaded_image_view
    void ReleaseData();

    // instance getters
    const MImageInfo& image_info() const
    { return image_info_; }

    MImageInfo& image_info()
    { return image_info_; }

    const MImageView& selected_view() const { return selected_view_; }

    const MImageView& loaded_view() const { return loaded_view_; }

    std::vector<int> loaded_dims() const override
    { return loaded_view_.dims(); }

    std::vector<int> loaded_xyz_dims() const override
    { return loaded_view_.xyz_dims(); }

    int loaded_n_channels() const override
    { return loaded_view_.n_channels(); }

    VoxelType loaded_voxel_type() const override
    { return loaded_view_.voxel_type(); }

private:
    /// allocate memory for selected view
    /// clears previous data and loaded_image_view_ by calling ClearData
    /// if the selected and loaded views have identical dimensions and voxel
    /// type, memory is not reallocated but loaded_image_view_ is still cleared
    void AllocateSelection() override;

    template <typename VType>
    bool LoadedDataIsType() const;

    uint8_t* data_impl(int volume_index) override;

    const uint8_t* const data_impl(int volume_index) const override;

    MImageInfo image_info_;
    /// these views are in the global context
    MImageView selected_view_, loaded_view_;

};

}

template <typename VType>
VType* mcp3d::MImage::Volume(int channel_number, int time)
{
    MCP3D_ASSERT(!empty())
    // during image read, its possible for memory to be allocated
    // based on selection, while loaded view is empty before read completes
    const mcp3d::MImageView& view = loaded_view_.empty() ? selected_view_ : loaded_view_;
    MCP3D_ASSERT(view.view_volume_index(channel_number) >= 0)
    MCP3D_ASSERT(time >= 0 && time < view.n_times())
    return reinterpret_cast<VType*>(data_[channel_number].get());
}

template <typename VType>
const VType* mcp3d::MImage::ConstVolume(int channel_number, int time) const
{
    MCP3D_ASSERT(!empty())
    // during image read, its possible for memory to be allocated
    // based on selection, while loaded view is empty before read completes
    const mcp3d::MImageView& view = loaded_view_.empty() ? selected_view_ : loaded_view_;
    MCP3D_ASSERT(view.view_volume_index(channel_number) >= 0)
    MCP3D_ASSERT(time >= 0 && time < view.n_times())
    int index = loaded_view_.view_volume_index(channel_number);
    return reinterpret_cast<const VType*>(data_[index].get());
}

template <typename VType>
VType* mcp3d::MImage::Plane(int channel_number, int view_z, int time)
{
    const mcp3d::MImageView& view = loaded_view_.empty() ? selected_view_ : loaded_view_;
    MCP3D_ASSERT(view.voxel_type() != mcp3d::VoxelType::UNKNOWN)
    return Volume<VType>(channel_number, time) + view.PlaneVoxels() * (long)view_z;
}

template <typename VType>
void mcp3d::MImage::AcquireData(std::vector<std::unique_ptr<VType[]>> &other_data)
{
    if (other_data.empty())
        return;
    data_.clear();
    if (!LoadedDataIsType<VType>())
        std::cout << "warning: currently loaded voxel type " << mcp3d::VoxelTypeToString(loaded_view_.voxel_type())
                  << ", acquiring voxel type " << mcp3d::VoxelTypeToString(mcp3d::TypeToVoxelType<VType>()) << std::endl;
    for (size_t i = 0; i < other_data.size(); ++i)
    {
        std::unique_ptr<uint8_t[]> other_ptr(reinterpret_cast<uint8*>(other_data[i].release()));
        data_.push_back(std::move(other_ptr));
    }
}

template <typename VType>
void mcp3d::MImage::AcquireData(std::unique_ptr<VType[]> &other_data, int volume_index)
{
    if (!other_data)
        return;
    if (volume_index != -1)
        MCP3D_ASSERT(volume_index >= 0 && volume_index < (int)(data_.size()))
    if (!LoadedDataIsType<VType>())
        std::cout << "warning: currently loaded voxel type " << mcp3d::VoxelTypeToString(loaded_view_.voxel_type())
                  << ", acquiring voxel type " << mcp3d::VoxelTypeToString(mcp3d::TypeToVoxelType<VType>()) << std::endl;
    std::unique_ptr<uint8_t[]> other_ptr(reinterpret_cast<uint8*>(other_data.release()));
    if (volume_index == -1)
        data_.push_back(std::move(other_ptr));
    else
        data_[volume_index] = std::move(other_ptr);
}

template <typename  VType>
bool mcp3d::MImage::LoadedDataIsType() const
{
    static_assert(std::is_arithmetic<VType>(), "can only pass arithmetic type template parameter");;
    return mcp3d::TypeToVoxelType<VType>() == loaded_view_.voxel_type();
}


#endif //MCP3D_MCP3D_IMAGE_HPP
