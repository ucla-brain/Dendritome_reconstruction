//
// Created by muyezhu on 12/1/18.
//

#ifndef MCP3D_MCP_TIFF_TILE_FORMAT_HPP
#define MCP3D_MCP_TIFF_TILE_FORMAT_HPP

#include "common/mcp3d_utility.hpp"
#include "mcp3d_image_formats.hpp"

namespace mcp3d
{

class MChannelPyrInfo;
class MChannelInfo;
class MImage;

/// this class supports uniform tiff image stack chunks that represent the
/// volume. does not enforce the xml tags of ometiff, but will name images as
/// .ome.tif to differentiate output from MTiffFormat
/// use *z[0-9]_y[0-9]_x[0-9].ome.tif to indicate the global coordinate of top
/// left voxel in first directory. each pyramid level will maintain identical
/// chunk xyz dimensions. padding is provided if image dimensions are not
/// wholely divisible by chunk dimensions during writing. reading operation
/// can not discern padding vs actual data
class MOmeTiffFormat: public MImageFormats
{
public:
    using json = nlohmann::json;

    MOmeTiffFormat(): MImageFormats(FileFormat::OMETIFF) {};

    bool CanRead() override { return true; }

    bool CanWrite() override { return true; }

    bool CanReadPartiallyComplete() override { return true; }

    bool CanWriteInChunk() override { return false; }

    bool CanWriteInChunkParallel() override { return false; }

    /// will sort img_paths
    MChannelPyrInfo ReadChannelPyrInfo(const MChannelPyrSlices& channel_pyr_slices, int resolution_level) override;

    void ValidateChannelInfo(const mcp3d::MChannelInfo &channel_info) override {}

    void WriteViewVolume(const MImage &image, const std::string &out_dir, const std::string &image_name_prefix) override {}

    void WriteImagePyramid(MImage &image, int channel, int parent_level, bool multi_threading) override
    { MCP3D_MESSAGE("not implemented") }

    #if MCP3D_MPI_BUILD

    void WriteImagePyramidMPI(MImage &image, int channel, int parent_level,
                              bool abort_all_on_fail,
                              const std::string &err_log_path,
                              const std::string &out_log_path,
                              MPI_Comm comm_writer) override {}
    #endif


private:
    std::string ImageNameSuffix(int z, int y, int x) const noexcept override
    { return "z" + PadNumStr(z, VOLUME_DIM_WIDTH) + "_y" + PadNumStr(y, VOLUME_DIM_WIDTH) + "_x" + PadNumStr(x, VOLUME_DIM_WIDTH) + ".ome.tif"; }

    std::vector<int> OptimalReadXyzDims(const MChannelPyrInfo &pyr_info) const override;

    void ReadChannelDataImpl(MImage &image, int channel_number) override;
};

}

#endif //MCP3D_MCP_TIFF_TILE_FORMAT_HPP
