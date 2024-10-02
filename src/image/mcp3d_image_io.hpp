//
// Created by muyezhu on 2/12/18.
//
#include <iostream>
#include <limits>
#include <memory>
#include <type_traits>
#include <map>
#include "mcp3d_image_info.hpp"
#include "mcp3d_image_view.hpp"
#include "mcp3d_image.hpp"
#include "mcp3d_image_formats.hpp"

#ifndef MCP3D_MCP3D_IMAGE_IO_HPP
#define MCP3D_MCP3D_IMAGE_IO_HPP

namespace mcp3d
{

class MImageIO
{
public:
    MImageIO();

    /// call MChannelLayout::RefreshChannelLayout and ReadChannelPyrInfos
    /// if pyr_level is negative, for all pyr levels found under channel_info.channel_layout_, read MChannelPyrInfo.
    /// otherwise only read MChannelPyrInfo at pyr_level
    /// call MChannelInfo::AssertValid
    /// pre exisitng MChannelPyrInfo for channel_name at pyr_level will be replaced by ReadChannelPyrInfos
    void ReadChannelInfo(MImage& image, const std::string& channel_name, int pyr_level = -1, bool load_from_json = true);

    void ReadChannelInfo(MImage& image, int channel_number, int pyr_level = -1, bool ignore_saved = false)
    { MCP3D_TRY(ReadChannelInfo(image, image.image_info().volume_layout().channel_name(channel_number, 0), pyr_level, ignore_saved);) }

    /// fill all volumes with zero before calling reading operations from specific formats
    /// this allows the portion of view outside of the selected level to be filled with black, as well as
    /// performing zero filling for partially written volume where view requests data that does not exist on disk yet
    /// if selected view is entirely out of volume bounds at selected level, nothing more will be done after zero filling
    void ReadData(MImage &image);

    void WriteViewVolume(const MImage &image, const std::string &out_dir, const std::string &img_name_prefix, FileFormat write_format);

    /// ui should ask user if over writing existing pyramid levels.
    /// this function will assume it should delete existing levels and make anew
    /// parent level and child level need not have the same ImageFormat
    /// at parent level the ImageFormats class must support reading, at child
    /// level the ImageFormats class must support writing
    void WriteImagePyramid(MImage &image, int channel, int parent_level, bool multi_threading, FileFormat write_format);

    bool ReadableFormat(FileFormat format);

    bool WritableFormat(FileFormat format);

    #if MCP3D_MPI_BUILD
    void WriteImagePyramidMPI(const std::string &img_root_dir, int channel,
                              int parent_level, FileFormat write_format,
                              bool abort_all_on_fail, MPI_Comm writer_comm);

    #endif

private:
    /// remove pre exisitng MChannelPyrInfo for channel_name at pyr_levels
    /// if load_from_json is true, call MChannelInfo::Load. this can lead to levels not included in pyr_levels to be loaded
    /// if MChannelPyrInfo at pyr_level has been loaded from json, and its format is same as the format found by
    /// channel_info.channel_layout_ for pyr_level, skip reading MChannelPyrInfo from disk image files at pyr_level
    void ReadChannelPyrInfos(MChannelInfo &channel_info, std::vector<int> pyr_levels, bool load_from_json);

    /// fill view entirely with zeros
    template<typename VType>
    static void FillViewWithZeros(MImage &image);

    std::map<FileFormat, std::unique_ptr<MImageFormats>> io_formats_;
};

}

template<typename VType>
void mcp3d::MImageIO::FillViewWithZeros(MImage &image)
{
    VType background_value = mcp3d::VoxelTypeZeroValue<VType>();
    for (int view_channel: image.selected_view().view_channels())
        mcp3d::SetConstantBlock<VType>(image.Volume<VType>(view_channel), image.selected_view().xyz_dims(), mcp3d::MImageBlock{}, background_value);
}


#endif //MCP3D_MCP3D_IMAGE_IO_HPP
