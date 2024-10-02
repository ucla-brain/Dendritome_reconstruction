//
// Created by muyezhu on 5/30/19.
//

#ifndef MCP3D_MCP3D_IMAGE_OPERATIONS_HPP
#define MCP3D_MCP3D_IMAGE_OPERATIONS_HPP

#include <type_traits>
#include <limits>
#include <memory>
#include <unsupported/Eigen/CXX11/Tensor>
#include "common/mcp3d_macros.hpp"
#include "common/mcp3d_types.hpp"
#include "common/mcp3d_utility.hpp"
#include "image_interface/mcp3d_voxel_types.hpp"
#include "image/mcp3d_image_utils.hpp"

namespace mcp3d
{

template <typename VType=uint16_t>
VType VolumeMin(void* data, const std::vector<int>& dims);

template <typename VType=uint16_t>
VType VolumeMax(void* data, const std::vector<int>& dims);

template <typename SrcVType = uint16_t, typename TargetVType = uint8_t>
void CastVolume(std::unique_ptr<uint8_t []>& src_data, const std::vector<int>& dims, std::unique_ptr<uint8_t []>& target_data);

// input is passed as uint8_t raw buffer. function should interpret the buffer
// as indicated by SrcVType and TargetVType
// floating point target types normalize to [0, 1]. integer types to [0, type_max]
template <typename SrcVType = uint16_t, typename TargetVType = uint8_t>
void NormalizeVolume(std::unique_ptr<uint8_t []>& src_data, const std::vector<int>& dims, std::unique_ptr<uint8_t []>& target_data);

// the q-th top percentile of an array of sorted values (min to max) is the
// q * array length-th value from back of the array, aka (1 - q) * array length-th
// value from front of the array. this index (1 - q) * array length is rounded to
// nearest integer. special: q = 0.0 -> max array value, q = 1.0 -> min array value
/// dims: data volume dimensions.
template<typename VType>
VType VolumeTopPercentile(void *data, const std::vector<int> &dims, double q);

/// block: MImageBlock within data volume
/// gives the top percentile value of the subvolume defined by block
template<typename VType>
VType VolumeBlockTopPercentile(void *data, const std::vector<int> &dims, const MImageBlock& block, double q);

/// data is asserted to be 3d volume. for the dimension selected by dim_id,
/// gives top percentile values of planes along the dimension
template <typename VType>
std::vector<VType> PlaneTopPercentiles(void *data, const std::vector<int> &dims, double q, int dim_id = 0);

/// map elements >= threshold_value to positive_value, and elements < threshold_value
/// to negative value
template <typename VType>
void ThresholdVolume(void* data, const std::vector<int> &dims, VType threshold_value, bool keep_values = false,
                     VType positive_value = std::numeric_limits<VType>::max(), VType negative_value = mcp3d::VoxelTypeZeroValue<VType>());


template <typename VType>
void ThresholdVolumeBlock(void* data, const std::vector<int> &dims, const MImageBlock& block, VType threshold_value, bool keep_values = false,
                          VType positive_value =std::numeric_limits<VType>::max(), VType negative_value = mcp3d::VoxelTypeZeroValue<VType>());

template <typename VType>
void ThresholdPlanes(void* data, const std::vector<int> &dims, const std::vector<VType>& threshold_values, int dim_id = 0,
                     bool keep_values = false, VType positive_value = std::numeric_limits<VType>::max(), VType negative_value = mcp3d::VoxelTypeZeroValue<VType>());

template <typename VType>
VType ThresholdVolumeTopPercentile(void* data, const std::vector<int> &dims, double q, bool keep_values = false,
                                   VType positive_value = std::numeric_limits<VType>::max(), VType negative_value = mcp3d::VoxelTypeZeroValue<VType>());

template <typename VType>
VType ThresholdVolumeBlockTopPercentile(void* data, const std::vector<int> &dims, const MImageBlock& block, double q, bool use_block_percentile = true,
                                        bool keep_values = false, VType positive_value = std::numeric_limits<VType>::max(),
                                        VType negative_value = mcp3d::VoxelTypeZeroValue<VType>());

template <typename VType>
std::vector<VType> ThresholdPlanesTopPercentile(void* data, const std::vector<int> &dims, double q, int dim_id = 0, bool keep_values = false,
                                                VType positive_value = std::numeric_limits<VType>::max(), VType negative_value = mcp3d::VoxelTypeZeroValue<VType>());

}

template <typename VType>
VType mcp3d::VolumeMin(void* data, const std::vector<int>& dims)
{
    static_assert(std::is_arithmetic<VType>::value, "VType must be arithmetic type");
    MCP3D_ASSERT(data)
    MCP3D_ASSERT(mcp3d::AllPositive(dims))
    mcp3d::MCPTensor1DMap<VType> src_map((VType*)data, mcp3d::ReduceProdSeq<size_t>(dims));
    mcp3d::MCPTensor0D<VType> src_min = src_map.template minimum();
    return src_min(0);
}

template <typename VType>
VType mcp3d::VolumeMax(void* data, const std::vector<int>& dims)
{
    static_assert(std::is_arithmetic<VType>::value, "VType must be arithmetic type");
    MCP3D_ASSERT(data)
    MCP3D_ASSERT(mcp3d::AllPositive(dims))
    mcp3d::MCPTensor1DMap<VType> src_map((VType*)data, mcp3d::ReduceProdSeq<size_t>(dims));
    mcp3d::MCPTensor0D<VType> src_max = src_map.template maximum();
    return src_max(0);
}

template <typename SrcVType, typename TargetVType>
void mcp3d::CastVolume(std::unique_ptr<uint8_t []>& src_data, const std::vector<int>& dims, std::unique_ptr<uint8_t []>& target_data)
{
    MCP3D_ASSERT(src_data)
    MCP3D_ASSERT(mcp3d::AllPositive(dims))
    bool in_place = &src_data == &target_data;
    bool same_type = std::is_same<SrcVType, TargetVType>::value;
    bool narrowing = (std::is_integral<SrcVType>() && std::is_integral<TargetVType>() && sizeof(SrcVType) > sizeof(TargetVType)) ||
                     (std::is_floating_point<SrcVType>() && std::is_integral<TargetVType>()) ||
                     (std::is_floating_point<SrcVType>() && std::is_floating_point<TargetVType>() && sizeof(SrcVType) > sizeof(TargetVType));
    if (narrowing)
    MCP3D_MESSAGE("warning: narrowing cast, may lose information")
    auto n = mcp3d::ReduceProdSeq<size_t>(dims);
    // if same type and cast in place, nothing needs to be done
    if (same_type && in_place)
        return;
    // if same type and not in place, copy
    else if (same_type)
    {
        target_data = std::make_unique<uint8_t []>(n * sizeof(TargetVType));
        memcpy(target_data.get(), src_data.get(), n * sizeof(TargetVType));
    }
    else
    {
        // make an intermediate buffer
        std::unique_ptr<uint8_t []> buffer(new uint8_t[n * sizeof(TargetVType)]);
        // make 1D tensormap to be dimension flexible
        mcp3d::MCPTensor1DMap<SrcVType> src_map((SrcVType*)(src_data.get()), n);
        // using tensor map instead of tensor here so eigen will not manage buffer memory
        mcp3d::MCPTensor1DMap<TargetVType> buffer_map((TargetVType*)(buffer.get()), n);
        buffer_map = src_map.template cast<TargetVType>();
        // if result expected in place, release buffer to src. otherwise release
        // buffer to target
        // operator = takes rvalue. same as src_data.reset(target_data.release())
        if (in_place)
            src_data = std::move(buffer);
        else
            target_data = std::move(buffer);
    }
};

template <typename SrcVType, typename TargetVType>
void mcp3d::NormalizeVolume(std::unique_ptr<uint8_t []>& src_data,
                            const std::vector<int>& dims,
                            std::unique_ptr<uint8_t []>& target_data)
{
    static_assert(std::is_arithmetic<SrcVType>() && std::is_arithmetic<TargetVType>(),
                  "must have arithmetic element types");
    MCP3D_ASSERT(src_data)
    MCP3D_ASSERT(mcp3d::AllPositive(dims))
    // find min and max
    auto src_min = (double)mcp3d::VolumeMin<SrcVType>(src_data.get(), dims),
         src_max = (double)mcp3d::VolumeMax<SrcVType>(src_data.get(), dims);
    // cast src to double buffer
    auto n = mcp3d::ReduceProdSeq<size_t>(dims);
    std::unique_ptr<uint8_t []> buffer;
    if (!std::is_same<SrcVType, double>::value)
        mcp3d::CastVolume<SrcVType, double>(src_data, dims, buffer);
    else
    {
        buffer = std::make_unique<uint8_t []>(n * sizeof(double));
        memcpy(buffer.get(), src_data.get(), n * sizeof(double));
    }

    // normalize range of buffer
    mcp3d::MCPTensor1DMap<double> buffer_map((double*)(buffer.get()), n);
    // min is 0
    buffer_map = buffer_map - src_min;
    src_max -= src_min;
    double multiplier, epsilon = 1e-12;
    if (std::is_floating_point<TargetVType>::value)
        multiplier = 1.0 / (src_max + epsilon);
    else
        multiplier = (double)(std::numeric_limits<TargetVType>::max()) / (src_max + epsilon);
    buffer_map  = buffer_map * multiplier;
    // round and check bounds
    buffer_map.round();
    buffer_map = buffer_map.cwiseMax(0.0);
    buffer_map = buffer_map.cwiseMin((double)(std::numeric_limits<TargetVType>::max()));
    // cast buffer to target type
    std::unique_ptr<uint8_t []> result;
    mcp3d::CastVolume<double, TargetVType>(buffer, dims, result);
    // release result to src if in place, otherwise to target
    bool in_place = &src_data == &target_data;
    if (in_place)
        src_data = std::move(result);
    else
        target_data = std::move(result);
};

template <typename VType>
VType mcp3d::VolumeTopPercentile(void *data, const std::vector<int> &dims, double q)
{
    static_assert(std::is_arithmetic<VType>::value, "VType must be arithmetic type");
    MCP3D_ASSERT(data)
    MCP3D_ASSERT(mcp3d::AllPositive(dims))
    auto n = mcp3d::ReduceProdSeq<long>(dims);
    MCP3D_ASSERT(q >= 0 && q <= 1)
    std::unique_ptr<VType []> data_copy(new (std::nothrow) VType[n]);
    MCP3D_ASSERT(data_copy);
    memcpy(data_copy.get(), data, sizeof(VType) * n);
    mcp3d::QuickSort<VType>(data_copy.get(), n);
    double index = (1 - q) * n;
    if (index < 0)
        return data_copy[0];
    if (index > n - 1)
        return data_copy[n - 1];
    auto low_index = (long)std::floor(index), high_index = (long)std::ceil(index);
    // nearest interpolation
    if (index - low_index <= high_index - index)
        return data_copy[low_index];
    else
        return data_copy[high_index];
}

template <typename VType>
VType mcp3d::VolumeBlockTopPercentile(void *data, const std::vector<int> &dims, const MImageBlock &block, double q)
{
    static_assert(std::is_arithmetic<VType>::value, "VType must be arithmetic type");
    MCP3D_ASSERT(dims.size() == 3 && mcp3d::AllPositive(dims))
    mcp3d::MImageBlock block_(block);
    for (int i = 0; i < 3; ++i)
        if (block_.extents_ptr()[i] == 0)
            block_.extents_ptr()[i] = dims[i] - block_.offsets_ptr()[i];
    if (block_.extents() == dims)
        return mcp3d::VolumeTopPercentile<VType>(data, dims, q);
    std::vector<int> strided_extents = mcp3d::StridedExtents(block_.extents(), block_.strides());
    auto n_elements = mcp3d::ReduceProdSeq<long>(strided_extents);
    std::unique_ptr<VType []> block_copy(new (std::nothrow) VType[n_elements]);
    mcp3d::CopyDataVolume<VType>(block_copy.get(), data, dims, block_);
    return mcp3d::VolumeTopPercentile<VType>(block_copy.get(), strided_extents, q);
}

template <typename VType>
std::vector<VType> mcp3d::PlaneTopPercentiles(void *data, const std::vector<int> &dims, double q, int dim_id)
{
    static_assert(std::is_arithmetic<VType>::value, "VType must be arithmetic type");
    MCP3D_ASSERT(dims.size() == 3)
    MCP3D_ASSERT(dim_id >= 0 && dim_id < 3)
    auto n_plane_elements = mcp3d::ReduceProdSeq<long>(dims) / (long)dims[dim_id];
    std::unique_ptr<VType []> plane_copy(new (std::nothrow) VType[n_plane_elements]);
    mcp3d::MImageBlock data_block({0, 0, 0}, {0, 0, 0}, {1, 1, 1});
    for (int i = 0; i < 3; ++i)
        data_block.extents_ptr()[i] = i == dim_id ? 1 : dims[i];
    std::vector<VType> plane_vals;
    for (int i = 0; i < dims[dim_id]; ++i)
    {
        data_block.offsets_ptr()[dim_id] = i;
        plane_vals.push_back(mcp3d::VolumeBlockTopPercentile<VType>(data, dims, data_block, q));
    }
    return plane_vals;
}

template <typename VType>
void mcp3d::ThresholdVolume(void *data, const std::vector<int> &dims, VType threshold_value, bool keep_values, VType positive_value, VType negative_value)
{
    static_assert(std::is_arithmetic<VType>::value, "VType must be arithmetic type");
    MCP3D_ASSERT(data)
    MCP3D_ASSERT(mcp3d::AllPositive(dims))
    auto n = mcp3d::ReduceProdSeq<long>(dims);
    auto vtype_data = (VType*)data;
    for (long i = 0; i < n; ++i)
        vtype_data[i] = vtype_data[i] >= threshold_value ?
                        (keep_values ? vtype_data[i] : positive_value) : negative_value;
}

template <typename VType>
void mcp3d::ThresholdVolumeBlock(void *data, const std::vector<int> &dims, const mcp3d::MImageBlock &block,
                                 VType threshold_value, bool keep_values, VType positive_value, VType negative_value)
{
    static_assert(std::is_arithmetic<VType>::value, "VType must be arithmetic type");
    MCP3D_ASSERT(data)
    MCP3D_ASSERT(dims.size() == 3 && mcp3d::AllPositive(dims))
    mcp3d::MImageBlock block_(block);
    for (int i = 0; i < 3; ++i)
        if (block_.extents_ptr()[i] == 0)
            block_.extents_ptr()[i] = dims[i] - block_.offsets_ptr()[i];
    if (block_.extents() == dims)
    {
        mcp3d::ThresholdVolume<VType>(data, dims, threshold_value, positive_value, negative_value);
        return;
    }
    auto vtype_data = (VType*)data;
    for (int z = block_.offsets_ptr()[0]; z < block_.offsets_ptr()[0] + block_.extents_ptr()[0]; z += block_.strides_ptr()[0])
    {
        for (int y = block_.offsets_ptr()[1]; y < block_.offsets_ptr()[1] + block_.extents_ptr()[1]; y += block_.strides_ptr()[1])
        {
            int64_t address = mcp3d::LinearAddressRobust(dims, z, y, block_.offsets_ptr()[2]);
            int x = block_.offsets_ptr()[2],
                xmax = block_.offsets_ptr()[2] + block_.extents_ptr()[2];
            while (x < xmax)
            {
                vtype_data[address] = vtype_data[address] >= threshold_value ?
                                      (keep_values ? vtype_data[address] : positive_value) : negative_value;
                address += block_.strides_ptr()[2];
                x += block_.strides_ptr()[2];
            }
        }
    }
}

template <typename VType>
void mcp3d::ThresholdPlanes(void *data, const std::vector<int> &dims, const std::vector<VType>& threshold_values, int dim_id,
                            bool keep_values, VType positive_value, VType negative_value)
{
    static_assert(std::is_arithmetic<VType>::value, "VType must be arithmetic type");
    MCP3D_ASSERT(dim_id >= 0 && dim_id < (int)dims.size())
    MCP3D_ASSERT(threshold_values.size() == (size_t)dims[dim_id])
    mcp3d::MImageBlock data_block({0, 0, 0}, {1, 1, 1}, {1, 1, 1});
    for (int i = 0; i < 3; ++i)
        data_block.extents_ptr()[i] = i == dim_id ? 1 : dims[i];
    for (int i = 0; i < dims[dim_id]; ++i)
    {
        data_block.offsets_ptr()[dim_id] = i;
        mcp3d::ThresholdVolumeBlock<VType>(data, dims, data_block, threshold_values[i], keep_values, positive_value, negative_value);
    }
}

template <typename VType>
VType mcp3d::ThresholdVolumeTopPercentile(void *data, const std::vector<int> &dims, double q, bool keep_values, VType positive_value, VType negative_value)
{
    static_assert(std::is_arithmetic<VType>::value, "VType must be arithmetic type");
    VType val = mcp3d::VolumeTopPercentile<VType>(data, dims, q);
    mcp3d::ThresholdVolume<VType>(data, dims, val, keep_values, positive_value, negative_value);
    return val;
}

template <typename VType>
VType mcp3d::ThresholdVolumeBlockTopPercentile(void *data, const std::vector<int> &dims, const MImageBlock &block, double q, bool use_block_percentile,
                                               bool keep_values, VType positive_value, VType negative_value)
{
    static_assert(std::is_arithmetic<VType>::value, "VType must be arithmetic type");
    mcp3d::MImageBlock block_(block);
    for (int i = 0; i < 3; ++i)
        if (block_.extents_ptr()[i] == 0)
            block_.extents_ptr()[i] = dims[i] - block_.offsets_ptr()[i];
    if (block_.extents() == dims)
        return mcp3d::ThresholdVolumeTopPercentile<VType>(data, dims, q, positive_value, negative_value);
    VType threshold_value;
    if (!use_block_percentile)
        threshold_value = mcp3d::VolumeTopPercentile<VType>(data, dims, q);
    else
        threshold_value = mcp3d::VolumeBlockTopPercentile<VType>(data, dims, block, q);
    mcp3d::ThresholdVolumeBlock<VType>(data, dims, block, threshold_value, keep_values, positive_value, negative_value);
    return threshold_value;
}

template <typename VType>
std::vector<VType> mcp3d::ThresholdPlanesTopPercentile(void *data, const std::vector<int> &dims, double q, int dim_id, bool keep_values,
                                                       VType positive_value, VType negative_value)
{
    static_assert(std::is_arithmetic<VType>::value, "VType must be arithmetic type");
    MCP3D_ASSERT(dim_id >= 0 && dim_id < (int)dims.size())
    std::vector<VType> threshold_values(mcp3d::PlaneTopPercentiles<VType>(data, dims, q, dim_id));
    mcp3d::ThresholdPlanes<VType>(data, dims, threshold_values, dim_id, keep_values, positive_value, negative_value);
    return threshold_values;
}

#endif //MCP3D_MCP3D_IMAGE_OPERATIONS_HPP
