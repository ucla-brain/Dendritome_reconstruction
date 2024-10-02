//
// Created by muyezhu on 3/13/19.
//

#ifndef MCP3D_MCP3D_HDF5_UTIL_HPP
#define MCP3D_MCP3D_HDF5_UTIL_HPP

#include <memory>
#include <vector>
#include <hdf5.h>
#include "mcp3d_voxel_types.hpp"

/// will close hid_t handles created by the function, but will not close handles received as argument
namespace mcp3d
{

hid_t Hdf5Handle(const std::string &hdf5_path);

std::vector<std::string> ObjectNamesInGroup(hid_t group_id);

void CloseHdfObject(hid_t object_id);

H5O_type_t HdfObjectType(hid_t object_id);

// retrieve attribute value into buffer. returns number of elements in
// attribute value. caller should know the attribute value element datatype
// in order to interpret content in buffer
int Hdf5AttributeValue(hid_t object_id, const char *attribute_name, std::unique_ptr<uint8_t[]> &buffer);

/// retrieve attribute value. interpret the buffer as uchar and return a corresponding string
std::string Hdf5AttributeString(hid_t object_id, const char *attribute_name);

bool IsHdf5Dataset(hid_t dataset_id);

bool IsChunkedHdf5Dataset(hid_t dataset_id);

VoxelType Hdf5DatasetVoxelType(hid_t dataset_id);

/// return 0 if dataset_id is not a dataset id
/// return 0 if dataspace is not simple
int Hdf5DatasetRank(hid_t dataset_id);

/// return empty vector if dataset_id is not a dataset id
/// return empty vector if dataspace is not simple
std::vector<int> Hdf5DatasetDimensions(hid_t dataset_id);

/// return empty vector if dataset_id is not a dataset id
/// dimensions of dataset chunk. if dataset is not chunked, return zero filled vector
std::vector<int> Hdf5DatasetChunkDimensions(hid_t dataset_id);

bool Hdf5ZlibFilterAvailable();

void SetHdf5DatasetZlibDeflate(hid_t dataset_id, int deflation);

hid_t DatasetCreationPropertyHandle(hid_t dataset_id);

hid_t MemoryDataSpace(hsize_t* mem_ds_dims, hsize_t* mem_ds_starts, hsize_t* mem_ds_counts, int rank = 3);

hid_t FileDataSpace(hid_t dataset_id, hsize_t* file_ds_starts, hsize_t* file_ds_strides, hsize_t* file_ds_counts);

/// will not free pointers
void ReadHdf5Dataset(hid_t dataset_id, hsize_t* file_ds_starts, hsize_t* file_ds_strides, hsize_t* ds_counts,
                     hsize_t* mem_ds_dims, hsize_t* mem_ds_starts, hid_t mem_type_id, void* buffer);

}

#endif //MCP3D_MCP3D_HDF5_UTIL_HPP
