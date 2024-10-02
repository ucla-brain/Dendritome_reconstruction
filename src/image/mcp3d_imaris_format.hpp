//
// Created by muyezhu on 1/30/19.
//

#ifndef MCP3D_MCP3D_HDF5_FORMAT_HPP_HPP
#define MCP3D_MCP3D_HDF5_FORMAT_HPP_HPP

#include <hdf5.h>
#include "common/mcp3d_common.hpp"
#include "image_interface/mcp3d_file_formats.hpp"
#include "mcp3d_image_formats.hpp"

/// the hdf5 format will support reading .ims, reading and writing .hdf5 files,
/// accessed by a single virtual hdf5 file. this allows parallel processing
/// and easier syntax. the virtual file should have extension .virtual.hdf5
/// note that the .ims files are often assumed to have a single time point 0
namespace mcp3d
{

class MImage;
class MChannelPyrSlices;
class MChannelPyrInfo;
class MChannelInfo;

class ImarisResolutionLevelInfo
{
public:
    ImarisResolutionLevelInfo(): imaris_path_(std::string{}), time_(-1), resolution_level_(-1), image_xyz_sizes_(std::vector<int>{}),
                                 chunk_xyz_dims_(std::vector<int>{}), voxel_type_(VoxelType::UNKNOWN) {}

    ImarisResolutionLevelInfo(const std::string &imaris_path, int resolution_level, int time = 0);

    int time() const
    { return time_; }

    int resolution_level() const
    { return resolution_level_; }

    bool is_chunked() const
    { return !chunk_xyz_dims_.empty(); }

    std::vector<int> image_xyz_sizes() const
    { return image_xyz_sizes_; }

    std::vector<int> chunk_xyz_dims() const
    { return chunk_xyz_dims_; }

    VoxelType voxel_type() const
    { return voxel_type_; }

private:
    std::string imaris_path_;
    int time_, resolution_level_;
    /// the chunk here is the hdf5 internal chunk, not file chunk size in MChannelInfo
    /// if the imaris file doesn't have chunked dataset layout, chunk_xyz_dims_ is empty
    /// note that image_xyz_sizes_[0] is the volume zdim
    std::vector<int> image_xyz_sizes_, chunk_xyz_dims_;
    VoxelType voxel_type_;
};

class MImarisFormat: public MImageFormats
{
public:
    MImarisFormat(): MImageFormats(FileFormat::IMARIS)  {};

    bool CanRead() override { return true; }

    bool CanWrite() override { return false; }

    /// imaris format does not support reading partially written volumes
    bool CanReadPartiallyComplete() override { return false; }

    bool CanWriteInChunk() override { return false; }

    bool CanWriteInChunkParallel() override { return false; }

    MChannelPyrInfo ReadChannelPyrInfo(const MChannelPyrSlices& channel_pyr_slices, int resolution_level) override;

    void ValidateChannelInfo(const mcp3d::MChannelInfo &image_info) override {};

    void WriteViewVolume(const MImage &img, const std::string &out_dir, const std::string &img_name_prefix) override {};

    /// imaris format should not be calling this function
    void WriteImagePyramid(MImage &image, int channel, int parent_level, bool multi_threading) override
    { MCP3D_RUNTIME_ERROR("not implemented") };

#if MCP3D_MPI_BUILD

    virtual void WriteImagePyramidMPI(MImage &image, int channel, int parent_level,
                                      bool abort_all_on_fail,
                                      const std::string &err_log_path,
                                      const std::string &out_log_path,
                                      MPI_Comm comm_writer) override {};

#endif

private:
    std::string ImageNameSuffix(int z, int y, int x) const noexcept override
    { return ".ims"; }

    /// if the imaris file is not chunked {zdim_, ydim_, xdim_} of pyr_info is returned.
    /// otherwise imaris chunk dimensions are returned
    std::vector<int> OptimalReadXyzDims(const MChannelPyrInfo &pyr_info) const override;

    void ReadChannelDataImpl(MImage &image, int channel_number) override;
};

}

#endif //MCP3D_MCP3D_HDF5_FORMAT_HPP_HPP
