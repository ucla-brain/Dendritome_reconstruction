//
// Created by muyezhu on 2/12/18.
//

#ifndef MCP3D_MCP3D_IMAGE_FORMATS_HPP
#define MCP3D_MCP3D_IMAGE_FORMATS_HPP

#include "common/mcp3d_common.hpp"
#include "mcp3d_channel_info.hpp"


namespace mcp3d
{

class MChannelPyrSlices;
class MChannelPyrInfo;
class MImage;

class MImageFormats
{
public:
    explicit MImageFormats(FileFormat format): format_(format) {}

    virtual ~MImageFormats() = default;

    virtual bool CanRead() = 0;

    virtual bool CanWrite() = 0;

    virtual bool CanReadPartiallyComplete() = 0;

    /// can modify data chunk without modifying data outside of the chunk
    virtual bool CanWriteInChunk() = 0;

    virtual bool CanWriteInChunkParallel() = 0;

    virtual MChannelPyrInfo ReadChannelPyrInfo(const MChannelPyrSlices& channel_pyr_slices, int resolution_level) = 0;

    /// this function reads MChannelPyrInfo for pyr_level volumes that are completely written
    /// some formats will check on the completeness assumption
    /// resolution_level argument is not required by tiff and ometiff formats (deducible from channel_pyr_dir) but
    /// needed for imaris format
    MChannelPyrInfo ReadChannelPyrInfo(const std::string &channel_pyr_dir, int resolution_level);

    /// assert MChannelPyrInfo at channel_name and pyr_level is present
    /// assert zyx is within the volume at channel_name level 0
    /// (z, y, x) is global coordinate (aka level 0). return image_path at pyr_level that contains (z, y, x)
    /// internally (z, y, x) will be mapped to (z', y', x') at pyr_level
    /// if image_path does not exist (e.g. in incomplete volumes), return empty string
    /// note that image_path does not exist is different from not having the required MChannelPyrInfo in
    /// image.image_info(), which causes the existence of image_path to be unknownable. it's also different
    /// from receiving (z, y, x) arguments that lies out side of level 0 volume, which is an invalid request
    std::string ImagePath(const MImage& image, const std::string& channel_name, int pyr_level, int z, int y = 0, int x = 0) const;

    std::string ImagePath(const MImage& image, int channel_number, int pyr_level, int z, int y = 0, int x = 0) const;

    /// assert image.image_info().HasChannelPyrInfo(channel_name, pyr_level)
    /// return the mininum zyx dimensions that can be read without waste io
    /// if selected view has offsets and extents both divisible by the optimal dimensions, reading is most efficient
    std::vector<int> OptimalReadXyzDims(const MImage& image, const std::string& channel_name, int pyr_level);

    /// validations common to all formats are performed by MImageIO.
    /// the validations performed by derived classes of MImageFormats are format
    /// specific
    virtual void ValidateChannelInfo(const mcp3d::MChannelInfo &image_info) = 0;

    /// view offsets must be in bounds of view level volume. views entirely out of bounds for the level
    /// should have been handled by MImageIO
    void ReadChannelData(MImage &image, int channel_number);

    virtual void WriteViewVolume(const MImage &image, const std::string &out_dir, const std::string &image_name_prefix) = 0;

    virtual void WriteImagePyramid(MImage &image, int channel_number, int parent_level, bool multi_threading) = 0;

    #if MCP3D_MPI_BUILD

    virtual void WriteImagePyramidMPI(MImage &image, int channel, int parent_level, bool abort_all_on_fail,
                                      const std::string &err_log_path, const std::string &out_log_path,
                                      MPI_Comm comm_writer) = 0;

    #endif

protected:
    /// return a string suffix of the image chunk with (z, y, x) and top left corner. (z, y, x) is local to whichever
    /// pyr level volume that requested this suffix. does not perform range checking for zyx
    virtual std::string ImageNameSuffix(int z, int y, int x) const noexcept = 0;

    virtual std::vector<int> OptimalReadXyzDims(const MChannelPyrInfo &pyr_info) const = 0;

    virtual void ReadChannelDataImpl(MImage &image, int channel_number) = 0;

    FileFormat format_;

private:
    /// as opposed to ImagePath, this function should understood zyx as local to pyr_level. (z, y, x) is assumed to have be asserted
    /// to have valid global correspondence.
    /// if (z, y, x) not in pur level volume and pyr_level > 0, return empty string
    /// if expected image path does not exist, return empty string
    std::string ImagePathIncompleteVolume(const MChannelPyrInfo &pyr_info, int z, int y, int x) const;

};

}

#if MCP3D_MPI_BUILD

#include <mpi.h>

#endif

#endif //MCP3D_MCP3D_IMAGE_FORMATS_HPP
