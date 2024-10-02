//
// Created by muyezhu on 1/30/19.
//

#ifndef MCP3D_IMS_IO_HPP_HPP
#define MCP3D_IMS_IO_HPP_HPP

#include <vector>
#include <hdf5.h>
#include "mcp3d_voxel_types.hpp"

namespace mcp3d
{
// will open and close file handle
// in general, functions close handles they opened within the function body,
// but not ones passes in as arguments
std::vector<std::string> ImarisResolutionNames(const std::string &imaris_path);

inline int NumberOfImarisResolutions(const std::string &imaris_path)
{ return (int)(ImarisResolutionNames(imaris_path).size()); }

// return sorted channel names
std::vector<std::string> ImarisResolutionChannelNames(hid_t imaris_id, int resolution_level, int time = 0);

// will open and close file handle
std::vector<std::string> ImarisResolutionChannelNames(const std::string &imaris_path, int resolution_level, int time = 0);

inline int ImarisResolutionNumberOfChannels(const std::string &imaris_path, int resolution_level, int time = 0)
{  return (int)(ImarisResolutionChannelNames(imaris_path, resolution_level, time).size()); }

// zyx dimensions of channel dataset without chunk padding
// will close file. the ImageSizeZ, ImageSizeY, ImageSizeX are represented
// as character arrays: e.g. '2', '0', '4', '8'
std::vector<int> ImarisChannelImageXyzSizes(hid_t imaris_id, int resolution_level, int channel_number, int time = 0);

// asserts image xyz sizes from all channels under the resolution level are equal
std::vector<int> ImarisResolutionImageXyzSizes(hid_t imaris_id, int resolution_level, int time = 0);

std::vector<int> ImarisChannelChunkXyzDims(hid_t imaris_id, int resolution_level, int channel_number, int time = 0);

std::vector<int> ImarisResolutionChunkXyzDims(hid_t imaris_id, int resolution_level, int time = 0);


VoxelType ImarisChannelVoxelType(hid_t imaris_id, int resolution_level, int channel_number, int time = 0);

VoxelType ImarisResolutionVoxelType(hid_t imaris_id, int resolution_level, int time = 0);

// handle to /DataSet/ResolutionLevel i
// will not close file
hid_t ImarisResolutionHandle(hid_t imaris_id, int resolution_level);

hid_t ImarisChannelHandle(hid_t imaris_id, int resolution_level, int channel_number, int time = 0);

hid_t ImarisTimePointHandle(hid_t imaris_id, int resolution_level, int time = 0);

// return handle to /DataSet/ResolutionLevel i/TimePoint j/Channel k/Data
hid_t ImarisDatasetHandle(hid_t imaris_id, int resolution_level, int channel_number, int time = 0);

}




#endif //MCP3D_IMS_IO_HPP_HPP
