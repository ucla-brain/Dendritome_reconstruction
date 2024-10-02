//
// Created by muyezhu on 7/19/19.
//

#include <exception>
#include "image_interface/mcp3d_file_formats.hpp"
#include "image_layout/mcp3d_channel_pyr_slices.hpp"
#include "image_layout/mcp3d_volume_layout.hpp"
#include "mcp3d_image.hpp"
#include "mcp3d_image_formats.hpp"

using namespace std;

mcp3d::MChannelPyrInfo mcp3d::MImageFormats::ReadChannelPyrInfo(const std::string &channel_pyr_dir, int resolution_level)
{
    MCP3D_TRY(return ReadChannelPyrInfo(mcp3d::MChannelPyrSlices{channel_pyr_dir}, resolution_level);)
}

string mcp3d::MImageFormats::ImagePath(const MImage &image, const string& channel_name, int pyr_level, int z, int y, int x) const
{
    try
    {
        if (!image.image_info().channel_pyr_info(channel_name, 0).ZyxInVolume(z, y, x))
            MCP3D_RUNTIME_ERROR("coordinates (" + to_string(z) + ", " + to_string(y) + ", " + to_string(x) + ") is out side of the global volume")
        if (!image.image_info().HasChannelPyrInfo(channel_name, pyr_level))
            MCP3D_RUNTIME_ERROR("no MChannelPyrInfo for channel " + channel_name + " and pyr level " + to_string(pyr_level) + " present in image.image_info()")
        vector<int> local_zyx = image.image_info().channel_info(channel_name).ZyxGlobalToLocal(pyr_level, z, y, x);
        const mcp3d::MChannelPyrInfo& pyr_info = image.image_info().channel_pyr_info(channel_name, pyr_level);
        if (pyr_info.volume_complete() == mcp3d::VolumeComplete::COMPLETE)
            return pyr_info.ImagePath(local_zyx[0], local_zyx[1], local_zyx[2]);
        return ImagePathIncompleteVolume(pyr_info, local_zyx[0], local_zyx[1], local_zyx[2]);
    }
    catch (...)
    {
        MCP3D_RETHROW(current_exception())
    }
}

string mcp3d::MImageFormats::ImagePath(const MImage &image, int channel_number, int pyr_level, int z, int y, int x) const
{
    MCP3D_TRY(return ImagePath(image, image.image_info().volume_layout().channel_name(channel_number, 0), pyr_level, z, y, x);)
}

vector<int> mcp3d::MImageFormats::OptimalReadXyzDims(const MImage &image, const std::string &channel_name, int pyr_level)
{
    if (!image.image_info().HasChannelPyrInfo(channel_name, pyr_level))
        MCP3D_RUNTIME_ERROR("no MChannelPyrInfo for channel " + channel_name + " and pyr level " + to_string(pyr_level) + " present in image.image_info()")
    return OptimalReadXyzDims(image.image_info().channel_pyr_info(channel_name, pyr_level));
}

void mcp3d::MImageFormats::ReadChannelData(MImage &image, int channel_number)
{
    try
    {
        if (image.selected_view().empty())
            MCP3D_RUNTIME_ERROR("can not read from image with empty view selection")
        if (image.image_info().channel_info(channel_number).file_format(image.selected_view().pyr_level()) != format_)
            MCP3D_RUNTIME_ERROR("file format of image at channel " + to_string(channel_number) + " pyr level " + to_string(image.selected_view().pyr_level()) +
                                " mismatches the instance")
        ReadChannelDataImpl(image, channel_number);
    }
    catch (...)
    {
        MCP3D_RETHROW(current_exception())
    }
}

string mcp3d::MImageFormats::ImagePathIncompleteVolume(const MChannelPyrInfo &pyr_info, int z, int y, int x) const
{
    MCP3D_ASSERT(pyr_info.volume_complete_ == mcp3d::VolumeComplete::INCOMPLETE);
    if (format_ == mcp3d::FileFormat::IMARIS)
        MCP3D_RUNTIME_ERROR("MImageFormats::ImagePathIncompleteVolume should not be called from IMARIS format")
    if (!pyr_info.ZyxInVolume(z, y, x))
        return string{};
    int zoffset = z - z % pyr_info.chunk_zdim_, yoffset = y - y % pyr_info.chunk_ydim_, xoffset = x - x % pyr_info.chunk_xdim_;
    string pyr_dir(pyr_info.channel_pyr_dir()),
            slice_name(pyr_info.channel_pyr_slices().SliceNameFromCoordinates(z, y, x)),
            image_name(ImageNameSuffix(zoffset, yoffset, xoffset));
    string image_path(mcp3d::JoinPath(pyr_dir, slice_name, image_name));
    return mcp3d::IsFile(image_path) ? image_path : string{};
}