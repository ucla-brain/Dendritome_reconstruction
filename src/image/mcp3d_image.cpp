//
// Created by muyezhu on 9/17/17.
//
#include <cstring>
#include <algorithm>
#include <opencv2/imgcodecs.hpp>
#include "common/mcp3d_utility.hpp"
#include "mcp3d_image_io.hpp"

using namespace std;

mcp3d::MImage::MImage(const MImage &other)
{

}

uint8_t* mcp3d::MImage::data_impl(int volume_index)
{
    if (volume_index < 0 or volume_index >= (int)data_.size())
        MCP3D_OUT_OF_RANGE("volume index out of range")
    return data_[volume_index].get();
}

const uint8_t* const mcp3d::MImage::data_impl(int volume_index) const
{
    if (volume_index < 0 or volume_index >= (int)data_.size())
        MCP3D_OUT_OF_RANGE("volume index out of range")
    return data_[volume_index].get();
}

void mcp3d::MImage::ReadImageInfo(const vector<int> &channel_numbers, int pyr_level, bool load_from_json)
{
    if (channel_numbers.empty())
    {
        MCP3D_MESSAGE("channel_numbers empty. no MChannelInfo will be read")
        return;
    }
    std::unique_ptr<MImageIO> io = make_unique<MImageIO>();
    try
    {
        for (int channel_number: channel_numbers)
            io->ReadChannelInfo(*this, channel_number, pyr_level, load_from_json);
    }
    catch (...)
    {
        MCP3D_PRINT_NESTED_EXCEPTION
        MCP3D_RETHROW(current_exception())
    }
}

void mcp3d::MImage::SelectView(const MImageBlock &view_block, const vector<int>& channel_numbers,
                               int pyr_level, bool interpret_block_as_local, bool allocate)
{
    vector<int> view_channels(channel_numbers);
    if (view_channels.empty())
        for (int i = 0; i < image_info_.volume_layout().n_channels(); ++i)
            view_channels.push_back(i);
    // read MChannelPyrInfo for view_channel at pyr_level if needed
    for (const auto& view_channel: view_channels)
    {
        if (!image_info_.HasChannelPyrInfo(view_channel, pyr_level))
        {
            ReadImageInfo(view_channel, pyr_level, true);
            if (pyr_level > 0 && !image_info().HasChannelPyrInfo(view_channel, 0))
                ReadImageInfo(view_channel, 0, true);
        }
    }
    // select view
    std::unique_ptr<MImageIO> io = make_unique<MImageIO>();
    try
    {
        selected_view_.SelectView(view_block, view_channels, pyr_level, interpret_block_as_local);
        if (allocate)
            AllocateSelection();
    }
    catch (...)
    {
        MCP3D_PRINT_NESTED_EXCEPTION
        MCP3D_RETHROW(current_exception())
    }
}

void mcp3d::MImage::ReadData(const string &mode)
{
    if (selected_view_.empty())
        MCP3D_RUNTIME_ERROR("no image view selected")
    if (selected_view_ == loaded_view_)
    {
        if (mode == "verbose")
            MCP3D_MESSAGE("selected image data is already loaded, do nothing")
        return;
    }
    if (mode == "verbose")
        selected_view_.PrintView();
    AllocateSelection();
    std::unique_ptr<MImageIO> io = make_unique<MImageIO>();
    try
    {
        io->ReadData(*this);
        loaded_view_.SelectView(selected_view_);
    }
    catch (...)
    {
        ClearLoaded();
        MCP3D_PRINT_NESTED_EXCEPTION
        MCP3D_RETHROW(current_exception())
    }
}

void mcp3d::MImage::AllocateSelection()
{
    if (selected_view_.empty())
        MCP3D_RUNTIME_ERROR("can only allocate memory for non empty image view selection")
    if (selected_view_.dims() == loaded_view_.dims() && selected_view_.voxel_type() == loaded_view_.voxel_type())
    {
        #ifdef VERBOSE
        MCP3D_MESSAGE("loaded image view have identical dimensions and voxel "
                      "type as selected view, no need to allocate memory")
        #endif
        return;
    }
    ClearLoaded();
    for (int i = 0; i < selected_view_.n_volumes(); ++i)
    {
        data_.push_back(make_unique<uint8_t[]>((size_t)selected_view_.VolumeBytes()));
        MCP3D_ASSERT(data_[i])
    }
}

void mcp3d::MImage::ClearLoaded()
{
    data_.clear();
    loaded_view_.Clear();
}

void mcp3d::MImage::SaveImageInfo()
{
    if (!mcp3d::MPIInitialized())
        // lambda function for macro expansion
        RANK0_CALL_SYNC([&](){image_info_.Save();})
    // in MPI execution consistent save should be controlled at higher level
    else
        image_info_.Save();
}

void mcp3d::MImage::WriteViewXYPlane(const std::string &img_path, int c, int z, int t)
{
    MCP3D_ASSERT(!data_.empty())
    MCP3D_ASSERT(z >= 0 && z < loaded_view_.zdim())
    MCP3D_ASSERT(c >= 0 && c < loaded_view_.n_channels())
    MCP3D_ASSERT(t >= 0 && t < loaded_view_.n_times())
    int cv_type = mcp3d::VoxelTypeToCVType(loaded_view_.voxel_type(), 1);
    uint8_t* ptr = Plane(c, z, t);
    cv::Mat m(loaded_view_.ydim(), loaded_view_.xdim(), cv_type, ptr);
    cv::imwrite(img_path, m);
}

void mcp3d::MImage::WriteViewVolume(const string &out_dir, const string& img_name_prefix, FileFormat volume_format)
{
    if (volume_format == FileFormat::UNKNOWN)
        MCP3D_DOMAIN_ERROR("view volume format can not be unknown")
    mcp3d::MakeDirectories(out_dir);
    if (loaded_view_.empty() || selected_view_.empty())
        cout << "no view is loaded or selected, do nothing" << endl;
    else if (loaded_view().empty())
        ReadData();

    unique_ptr<MImageIO> image_io = make_unique<MImageIO>();
    MCP3D_TRY(image_io->WriteViewVolume(*this, out_dir, img_name_prefix, volume_format);)
}

void mcp3d::MImage::WriteImagePyramids(int channel_number, int start_parent_level, int end_parent_level,
                                       bool multi_threading, bool save_image_info, FileFormat write_format)
{
    MCP3D_ASSERT(!image_info_.empty())
    if (!image_info_.volume_layout().HasChannel(channel_number))
        MCP3D_RUNTIME_ERROR("channel number " + to_string(channel_number) + " is not found under the MImage instance (volume_root_dir = " +
                            image_info_.volume_layout().volume_root_dir() + ")")
    if (!image_info_.volume_layout().HasChannelPyrLevel(channel_number, start_parent_level))
        MCP3D_RUNTIME_ERROR("pyr level " + to_string(start_parent_level) + " directory is not found under channel " + to_string(channel_number) +
                            " of the MImage instance (volume_root_dir = " + image_info_.volume_layout().volume_root_dir() + ")")
    if (!image_info_.HasChannelPyrInfo(channel_number, start_parent_level))
    {
        ReadImageInfo(channel_number, start_parent_level, false);
        if (!image_info_.HasChannelPyrInfo(channel_number, start_parent_level))
            MCP3D_RUNTIME_ERROR("pyr level " + to_string(start_parent_level) + " MChannelPyrInfo is not found under channel " + to_string(channel_number) +
                                " of the MImage instance (volume_root_dir = " + image_info_.volume_layout().volume_root_dir() + ")")
    }

    int mappable_level = image_info_.MaxMappablePyrLevel();
    end_parent_level =  end_parent_level > start_parent_level ? min(mappable_level, end_parent_level) : mappable_level;
    for (int i = start_parent_level; i < end_parent_level; ++i)
        MCP3D_TRY(WriteImagePyramid(channel_number, i, multi_threading, save_image_info, write_format);)
}

void mcp3d::MImage::WriteImagePyramid(int channel_number, int parent_level, bool multi_threading, bool save_image_info, FileFormat write_format)
{
    try
    {
        mcp3d::MImageIO image_io {};
        image_io.WriteImagePyramid(*this, channel_number, parent_level, multi_threading, write_format);
        image_io.ReadChannelInfo(*this, channel_number, parent_level + 1, false);
        if (save_image_info)
            SaveImageInfo();
    }
    catch (...)
    {
        MCP3D_RETHROW(current_exception())
    }
}

void mcp3d::MImage::ReleaseData()
{
    for (size_t i = 0; i < data_.size(); ++i)
        data_[i].release();
    data_.clear();
    loaded_view_.Clear();
}

