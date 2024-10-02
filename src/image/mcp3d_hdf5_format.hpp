//
// Created by muyezhu on 9/13/20.
//

#ifndef MCP3D_MCP3D_HDF5_FORMAT_HPP
#define MCP3D_MCP3D_HDF5_FORMAT_HPP

#include "common/mcp3d_common.hpp"
#include "image_interface/mcp3d_voxel_types.hpp"
#include "image_interface/mcp3d_file_formats.hpp"
#include "mcp3d_image_formats.hpp"

namespace mcp3d
{

class MChannelPyrSlices;
class MChannelPyrInfo;
class MChannelInfo;
class MImage;

class MHdf5Info
{
public:
    MHdf5Info(): hdf5_file_path_(std::string{}), xyz_dims_(std::vector<int>{}), chunk_xyz_dims_(std::vector<int>{}),
                 volume_xyz_dims_(std::vector<int>{}),voxel_type_(VoxelType::UNKNOWN) {}

    explicit MHdf5Info(const std::string& hdf5_file_path);

    std::vector<int> xyz_dims() const
    { return xyz_dims_; }

    std::vector<int> chunk_xyz_dims() const
    { return chunk_xyz_dims_; }

    std::vector<int> volume_xyz_dims() const
    { return volume_xyz_dims_; }

    VoxelType voxel_type() const
    { return voxel_type_; }

    bool is_chunked() const
    { return !chunk_xyz_dims_.empty(); }

private:
    std::string hdf5_file_path_;
    /// chunk dimensions refer to the HDF5 chunked storage concept.
    /// volume xyz dims is the dimension of an entire 3d volume, consisting of potentially multiple hdf5 files.
    /// all hdf5 subvolumes are of identical dimensions, therefore summing their total dimensions does not
    /// yiled volume dimensions. along each axis naxis subvolumes are created such that their total dimensions
    /// >= volume dimensions, and naxis is the smallest number that satisfy this requirement.
    /// volume dimensions are encoded as dataset properties
    std::vector<int> xyz_dims_, chunk_xyz_dims_, volume_xyz_dims_;
    VoxelType voxel_type_;
};

/// similar to ome tiff design. hdf5 file chunks with zyx chunk dimensions should
/// include volume zyx dimension in metadata fields
/// /"mcp3d_data": hdf5 dataset
///         attributes (attached to mcp3d_data): "volume xdim", "volume ydim", "volume zdim"
///                                              similar to imaris ImageSizeX
class MHdf5Format: public MImageFormats
{
public:
    MHdf5Format(): MImageFormats(FileFormat::HDF5) {};

    bool CanRead() override { return true; }

    bool CanWrite() override { return false; }

    bool CanReadPartiallyComplete() override { return true; }

    bool CanWriteInChunk() override { return true; }

    bool CanWriteInChunkParallel() override { return true; }

    MChannelPyrInfo ReadChannelPyrInfo(const MChannelPyrSlices& channel_pyr_slices, int resolution_level) override;

    void ValidateChannelInfo(const mcp3d::MChannelInfo &image_info) override {};

    void WriteViewVolume(const MImage &img, const std::string &out_dir, const std::string &img_name_prefix) override
    { MCP3D_MESSAGE("not implemented") };

    void WriteImagePyramid(MImage &image, int channel, int parent_level, bool multi_threading) override
    { MCP3D_MESSAGE("not implemented") };

#if MCP3D_MPI_BUILD

    virtual void WriteImagePyramidMPI(MImage &image, int channel, int parent_level,
                                      bool abort_all_on_fail,
                                      const std::string &err_log_path,
                                      const std::string &out_log_path,
                                      MPI_Comm comm_writer) override {};

#endif

    friend class MHdf5Info;

private:
    std::string ImageNameSuffix(int z, int y, int x) const noexcept override
    { return "z" + PadNumStr(z, VOLUME_DIM_WIDTH) + "_y" + PadNumStr(y, VOLUME_DIM_WIDTH) + "_x" + PadNumStr(x, VOLUME_DIM_WIDTH) + ".hdf5"; }

    /// if the hdf5 files are not chunked storage, {chunk_zdim_, chunk_ydim_, chunk_xdim_} of pyr_info is returned.
    /// otherwise hdf5 dataset chunk dimensions are returned
    std::vector<int> OptimalReadXyzDims(const MChannelPyrInfo &pyr_info) const override;

    /// returns dataset id for /dataset_name()
    static hid_t DatasetId(hid_t hdf5_id);

    static std::string dataset_name()
    { return "mcp3d_data"; }

    /// deos not assert axis to be known
    static std::string VolumeDimAttributeName(ChannelAxis axis)
    { return "volume_" + mcp3d::ChannelAxisStr(axis) + "dim"; }

    void ReadChannelDataImpl(MImage& image, int channel_number) override;

};

}

#endif //MCP3D_MCP3D_HDF5_FORMAT_HPP
