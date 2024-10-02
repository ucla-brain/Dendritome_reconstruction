//
// Created by muyezhu on 9/13/20.
//
#include <cstring>
#include <vector>
#include <hdf5.h>
#include "image_interface/mcp3d_hdf5_utils.hpp"
#include "image_layout/mcp3d_channel_pyr_slices.hpp"
#include "image_layout/mcp3d_channel_layout.hpp"
#include "mcp3d_channel_info.hpp"
#include "mcp3d_image.hpp"
#include "mcp3d_hdf5_format.hpp"

using namespace std;

mcp3d::MHdf5Info::MHdf5Info(const string& hdf5_file_path): hdf5_file_path_(hdf5_file_path)
{
    hid_t hdf5_id = mcp3d::Hdf5Handle(hdf5_file_path_);
    hid_t  dataset_id = mcp3d::MHdf5Format::DatasetId(hdf5_id);
    xyz_dims_ = mcp3d::Hdf5DatasetDimensions(dataset_id);
    if (xyz_dims_.size() != 3)
        MCP3D_RUNTIME_ERROR("expecting dataset wtih 3 dimensions, got " + to_string(xyz_dims_.size()) + " dimensions instead")
    chunk_xyz_dims_ = mcp3d::Hdf5DatasetChunkDimensions(dataset_id);
    if (chunk_xyz_dims_.size() != 3 && !chunk_xyz_dims_.empty())
        MCP3D_RUNTIME_ERROR("expecting dataset stored as rank 3 chunks or no chunk")
    voxel_type_ = mcp3d::Hdf5DatasetVoxelType(dataset_id);
    MCP3D_ASSERT(voxel_type_ != mcp3d::VoxelType::UNKNOWN)
    for (const auto axis: vector<mcp3d::ChannelAxis>({mcp3d::ChannelAxis::Z, mcp3d::ChannelAxis::Y, mcp3d::ChannelAxis::X}))
        volume_xyz_dims_.push_back(stoi(mcp3d::Hdf5AttributeString(dataset_id, MHdf5Format::VolumeDimAttributeName(axis).c_str())));
    H5Dclose(dataset_id);
    H5Fclose(hdf5_id);
}

mcp3d::MChannelPyrInfo mcp3d::MHdf5Format::ReadChannelPyrInfo(const mcp3d::MChannelPyrSlices &channel_pyr_slices, int resolution_level)
{
    // assert MChannelPyrSlices has HDF5 format
    MCP3D_ASSERT(channel_pyr_slices.file_format() == format_)
    // construct MChannelPyrInfo with channel_pyr_slices
    mcp3d::MChannelPyrInfo channel_pyr_info{channel_pyr_slices};
    MCP3D_ASSERT(mcp3d::MChannelLayout::DirPyrLevel(channel_pyr_info.channel_pyr_dir()) == resolution_level)
    channel_pyr_info.resolution_level_ = resolution_level;
    int n = channel_pyr_info.TotalNumberOfImageChunks();
    // because channel_pyr_slices.file_format() was determined to be HDF5, n should be positive
    MCP3D_ASSERT(n > 0)
    for (const auto& item: channel_pyr_info.slice_image_names_)
    {
        if (!item.second.empty())
        {
            mcp3d::MHdf5Info hdf5_info(item.second[0]);
            // dimensions of the volume at pyr level
            vector<int> xyz_dims(hdf5_info.volume_xyz_dims());
            channel_pyr_info.zdim_ = xyz_dims[0];
            channel_pyr_info.ydim_ = xyz_dims[1];
            channel_pyr_info.xdim_ = xyz_dims[2];
            /// dimensions of individual file chunk
            xyz_dims = hdf5_info.xyz_dims();
            channel_pyr_info.chunk_zdim_ = xyz_dims[0];
            channel_pyr_info.chunk_ydim_ = xyz_dims[1];
            channel_pyr_info.chunk_xdim_ = xyz_dims[2];
            channel_pyr_info.voxel_type_ = hdf5_info.voxel_type();
            break;
        }
    }
    return channel_pyr_info;
}

vector<int> mcp3d::MHdf5Format::OptimalReadXyzDims(const mcp3d::MChannelPyrInfo &pyr_info) const
{
    MCP3D_ASSERT(pyr_info.file_format_ == mcp3d::FileFormat::HDF5)
    for (const auto& item: pyr_info.slice_image_names_)
    {
        if (!item.second.empty())
        {
            mcp3d::MHdf5Info hdf5_info(item.second[0]);
            if (hdf5_info.is_chunked())
                return hdf5_info.chunk_xyz_dims();
            else
                return hdf5_info.xyz_dims();
        }
    }
    MCP3D_RUNTIME_ERROR("no files found in MChannelPyrInfo::slice_image_names_")
}

hid_t mcp3d::MHdf5Format::DatasetId(hid_t hdf5_id)
{
    hid_t  dataset_id = H5Dopen(hdf5_id, mcp3d::MHdf5Format::dataset_name().c_str(), H5P_DEFAULT);
    if (dataset_id < 0)
    MCP3D_RUNTIME_ERROR("can not open dataset \"" + mcp3d::MHdf5Format::dataset_name() + "\"")
    return dataset_id;
}

void mcp3d::MHdf5Format::ReadChannelDataImpl(MImage &image, int channel_number)
{
    const MImageView& view = image.selected_view();
    int pyr_level = image.selected_view().pyr_level();
    const MChannelPyrInfo& pyr_info = image.image_info().channel_pyr_info(channel_number, pyr_level);

    // determine maximum zyx positions in pyr level volume defined by selected view
    vector<int> xyz_maxs = mcp3d::Minimum(pyr_info.xyz_dims(), view.pyr_level_offsets() + view.pyr_level_extents());
    int xmax = xyz_maxs[2], ymax = xyz_maxs[1], zmax = xyz_maxs[0];
    // pyr level volume offsets, strides
    int volume_xoffset = view.pyr_level_offsets()[2], volume_yoffset = view.pyr_level_offsets()[1], volume_zoffset = view.pyr_level_offsets()[0];
    int volume_xstride = view.pyr_level_strides()[2], volume_ystride = view.pyr_level_strides()[1], volume_zstride = view.pyr_level_strides()[0];
    // hdf5 file dimensions
    int chunk_xdim = pyr_info.xyz_chunk_dims()[2], chunk_ydim = pyr_info.xyz_chunk_dims()[1], chunk_zdim = pyr_info.xyz_chunk_dims()[0];

    // traverse yx locations. at each location, read all z planes requested and available in file chunk
    int xcur = volume_xoffset, ycur, zcur;
    // tracks offset of the read opration into the view buffer
    int view_xoffset = 0, view_yoffset = 0, view_zoffset = 0;
    // tracks number of voxels to transfer to view buffer
    int view_xcount, view_ycount, view_zcount;
    // arguments to hdf5 library
    hsize_t file_ds_starts[3];  // zyxcur
    hsize_t file_ds_strides[3] = {(hsize_t)volume_zstride, (hsize_t)volume_ystride, (hsize_t)volume_xstride};
    hsize_t ds_counts[3];  // shared by mem and file space. view zyx counts
    hsize_t mem_ds_dims[3] = {(hsize_t)view.zdim(), (hsize_t)view.ydim(), (hsize_t)view.xdim()};
    hsize_t mem_ds_starts[3];  //view zyx offsets
    hid_t mem_type_id = mcp3d::VoxelTypeToH5Type(view.voxel_type());

    string hdf5_path;
    while (xcur < xmax)
    {
        ycur = volume_yoffset;
        view_yoffset = 0;
        int residual_stride;
        // number of voxels along x axis to be read from current file chunk, accounting for strides
        view_xcount = mcp3d::StridedExtent(min(chunk_xdim - xcur % chunk_xdim, xmax - xcur), volume_xstride);
        while (ycur < ymax)
        {
            zcur = volume_zoffset;
            view_zoffset = 0;
            view_ycount = mcp3d::StridedExtent(min(chunk_ydim - ycur % chunk_ydim, ymax - ycur), volume_ystride);
            while (zcur < zmax)
            {
                view_zcount = mcp3d::StridedExtent(min(chunk_zdim - zcur % chunk_zdim, zmax - zcur), volume_zstride);
                hdf5_path = ImagePath(image, channel_number, pyr_level, zcur, ycur, xcur);
                if (!hdf5_path.empty())
                {
                    hid_t hdf5_id = mcp3d::Hdf5Handle(hdf5_path);
                    hid_t dataset_id = mcp3d::MHdf5Format::DatasetId(hdf5_id);
                    file_ds_starts[0] = (hsize_t)zcur;
                    file_ds_starts[1] = (hsize_t)ycur;
                    file_ds_starts[2] = (hsize_t)xcur;
                    mem_ds_starts[0] = (hsize_t)view_zoffset;
                    mem_ds_starts[1] = (hsize_t)view_yoffset;
                    mem_ds_starts[2] = (hsize_t)view_xoffset;
                    ds_counts[0] = (hsize_t)view_zcount;
                    ds_counts[1] = (hsize_t)view_ycount;
                    ds_counts[2] = (hsize_t)view_xcount;
                    mcp3d::ReadHdf5Dataset(dataset_id, file_ds_starts, file_ds_strides, ds_counts,
                                           mem_ds_dims, mem_ds_starts, mem_type_id, image.Volume<uint8_t>(channel_number));
                    H5Dclose(dataset_id);
                    H5Fclose(hdf5_id);
                }
                residual_stride = mcp3d::ResidualStride(chunk_zdim - zcur % chunk_zdim, volume_zstride);
                // bring zcur to beginning for next file chunk, add residual stride
                zcur += chunk_zdim - zcur % chunk_zdim + residual_stride;
                zcur += 0; // todo
                view_zoffset += view_zcount;
            }
            residual_stride = mcp3d::ResidualStride(chunk_ydim - ycur % chunk_ydim, volume_ystride);
            ycur += chunk_ydim - ycur % chunk_ydim + residual_stride;
            view_yoffset += view_ycount;
        }
        residual_stride = mcp3d::ResidualStride(chunk_xdim - xcur % chunk_xdim, volume_xstride);
        xcur += chunk_xdim - xcur % chunk_xdim + residual_stride;
        view_xoffset += view_xcount;
    }
}


