//
// Created by muyezhu on 2/12/18.
//

#ifndef MCP3D_MCP3D_TIFF_FORMAT_HPP
#define MCP3D_MCP3D_TIFF_FORMAT_HPP

#include <opencv2/core/core.hpp>
#include <opencv2/imgcodecs/imgcodecs.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include "common/mcp3d_common.hpp"
#include "image_interface/mcp3d_file_formats.hpp"
#include "mcp3d_image_constants.hpp"
#include "mcp3d_image_macros.hpp"
#include "image_interface/mcp3d_tiff_io.hpp"
#include "mcp3d_image.hpp"
#include "mcp3d_image_formats.hpp"

namespace mcp3d
{

class MChannelPyrSlices;

/// this class supports sequence of tiff images of single channel and time point
/// aka each tiff image in the sequence is presumed to have a single directory
/// containing all data of one z plane
/// single image file that contain multiple z level, channels or time points
/// should be organized as tif stack or hdf
class MTiffFormat: public MImageFormats
{
public:
    MTiffFormat(): MImageFormats(FileFormat::TIFF)  {};

    bool CanRead() override { return true; }

    bool CanWrite() override { return true; }

    bool CanReadPartiallyComplete() override { return true; }

    bool CanWriteInChunk() override { return false; }

    bool CanWriteInChunkParallel() override { return false; }

    /// read a single pyramid level. the channel folder slices can only be non flat for z axis
    MChannelPyrInfo ReadChannelPyrInfo(const MChannelPyrSlices& channel_pyr_slices, int resolution_level) override;

    void ValidateChannelInfo(const MChannelInfo &channel_info) override;

    // does not allocation memory
    template <typename VType>
    void ReadPartialTiff(const std::string &tiff_path, VType* buffer, const MImage &image);

    void WriteViewVolume(const MImage &img, const std::string &out_dir, const std::string &img_name_prefix) override {};

    // divide image chunks to write evenly among threads. output is tiled tiff
    // maintain in memory tile strips of paren data
    // even if on disk storage of tiff image has rgb pixel types the pyramid will
    // be generated in grey
    void WriteImagePyramid(MImage &image, int channel_number, int parent_level, bool multi_threading) override;

    #if MCP3D_MPI_BUILD

    void WriteImagePyramidMPI(MImage &image, int channel_number, int parent_level,
                              bool abort_all_on_fail,
                              const std::string &err_log_path,
                              const std::string &out_log_path,
                              MPI_Comm comm_writer) override;

    #endif

private:
    /// tiff images does not need to include yx
    /// return z[0-9]{6}.tif
    std::string ImageNameSuffix(int z, int y, int x) const noexcept override
    { return "z" + PadNumStr(z, VOLUME_DIM_WIDTH) + ".tif"; }

    std::vector<int> OptimalReadXyzDims(const MChannelPyrInfo &pyr_info) const override;

    void ReadChannelDataImpl(MImage &image, int channel_number) override
    { MCP3D_TRY(IMAGE_SELECTED_TYPED_CALL(ReadChannelTiffDataImpl, image, image, channel_number))}

    template <typename VType>
    void ReadChannelTiffDataImpl(MImage &image, int channel_number);

    /// assert HasChannelPyrLevel amd HasChannelPyrInfo
    void WriteImagePyramidImpl(const mcp3d::MImage &image, int channel_number, int parent_level, int z_level);

};

}

// this method does not manage memory of MImage
template <typename VType>
void mcp3d::MTiffFormat::ReadChannelTiffDataImpl(mcp3d::MImage &image, int channel_number)
{
    const mcp3d::MChannelInfo& channel_info = image.image_info().channel_info(channel_number, 0);
    const mcp3d::MImageView& view = image.selected_view();

    // dimension of selected image view at current pyr_level, accounting for strides
    int view_xdim = view.xdim(), view_ydim = view.ydim();
    // dimension of selected image view at current pyr_level, assuming unit strides
    int view_xextent = view.pyr_level_extents()[2],
        view_yextent = view.pyr_level_extents()[1],
        view_zextent = view.pyr_level_extents()[0];
    int view_xstride = view.pyr_level_strides()[2],
        view_ystride = view.pyr_level_strides()[1],
        view_zstride = view.pyr_level_strides()[0];
    // the image portion covered by view assuming unit strides will be read first
    // if strides are not all unit, a transfer will happen according to the strides
    long n_plane_voxels = view.PlaneVoxels();
    std::unique_ptr<VType[]> ptr_temp;
    bool rgb_unit_stride = false;
    // used for intermediate buffer to only fill background once (if needed)
    bool background_filled = false;
    for (int z_image = view.pyr_level_offsets()[0], z_view = 0;
         z_image < view.pyr_level_offsets()[0] + view_zextent;
         z_image += view_zstride, ++z_view)
    {
        // if z value out of boundary, exist loop. background pixels are filled in MImageIO
        if (z_image >= channel_info.zdim(image.selected_view().pyr_level()))
            break;
        VType* ptr_img = image.Plane<VType>(channel_number, z_view);
        std::string img_path(ImagePath(image, channel_number, image.selected_view().pyr_level(), z_image, 0, 0));
        mcp3d::TiffDirectoryInfo tiff_info(img_path);
        // if unit stride and 1 sample per pixel in tiff image, read into ptr_img
        if (view.is_unit_strided("xy") && tiff_info.samples_per_pixel == 1)
        {
            ReadPartialTiff<VType>(img_path, ptr_img, image);
        }
        // otherwise read to an intermediate buffer, then transfer
        else
        {
            bool fill_background = image.selected_view().PartiallyOutOfPyrImageBoundary() && (!background_filled);
            // check if rgb, if so convert to gray before copy
            if (tiff_info.samples_per_pixel > 1)
            {
                long ptr_temp_addr = z_view * n_plane_voxels;
                if (z_image == view.pyr_level_offsets()[0])
                    ptr_temp = std::make_unique<VType[]>((size_t) view.VolumeVoxels());
                if (fill_background)
                {
                    background_filled = true;
                    memset(ptr_temp.get(), 0, (size_t) view.VolumeVoxels() * (size_t) view.VoxelBytes());
                }
                MCP3D_ASSERT(ptr_temp.get())
                std::unique_ptr<VType[]> ptr_cv = std::make_unique<VType[]>((size_t)view_xextent * (size_t)view_yextent * (size_t)tiff_info.samples_per_pixel);
                MCP3D_ASSERT(ptr_cv.get())
                ReadPartialTiff<VType>(img_path, ptr_cv.get(), image);
                cv::Mat m_in(view_yextent, view_xextent, mcp3d::TypeToCVType<VType>(tiff_info.samples_per_pixel), ptr_cv.get());
                cv::Mat m_out(view_yextent, view_xextent, mcp3d::TypeToCVType<VType>(1), ptr_temp.get() + ptr_temp_addr);
                SERIALIZE_OPENCV_MPI
                #if CV_MAJOR_VERSION < 4
                    cv::cvtColor(m_in, m_out, CV_RGB2GRAY);
                #else
                    cv::cvtColor(m_in, m_out, cv::COLOR_RGB2GRAY);
                #endif
                // if unit stride, do unique pointer ownership transfer
                if (view.is_unit_strided("xy"))
                    rgb_unit_stride = true;
            }
            else
            {
                size_t n_elements = (size_t)view_xextent * (size_t)view_yextent * (size_t)tiff_info.samples_per_pixel;
                ptr_temp = std::make_unique<VType[]>(n_elements);
                MCP3D_ASSERT(ptr_temp.get())
                if (fill_background)
                {
                    background_filled = true;
                    memset(ptr_temp.get(), 0, n_elements * view.VoxelBytes());
                }
                ReadPartialTiff<VType>(img_path, ptr_temp.get(), image);
                mcp3d::CopyDataVolume<VType>(ptr_img, {1, view_ydim, view_xdim}, ptr_temp.get(),
                                             {1, view_yextent, view_xextent}, mcp3d::MImageBlock{},
                                             mcp3d::MImageBlock({}, {}, {1, view_ystride, view_xstride}));
            }
        }
    }
    // unique pointer ownership transfer
    if (rgb_unit_stride)
        image.AcquireData<VType>(ptr_temp, channel_number);
}

template <typename VType>
void mcp3d::MTiffFormat::ReadPartialTiff(const std::string &tiff_path, VType *buffer, const mcp3d::MImage &image)
{
    MCP3D_ASSERT(buffer)
    if (!mcp3d::IsFile(tiff_path))
        MCP3D_OS_ERROR(tiff_path + " is not a file")
    const mcp3d::MImageView& view = image.selected_view();
    mcp3d::TransferToSubimg(tiff_path, buffer, view.pyr_level_extents()[1], view.pyr_level_extents()[2], view.pyr_level_offsets()[2], view.pyr_level_offsets()[1]);
}



#endif //MCP3D_MCP3D_TIFF_FORMAT_HPP
