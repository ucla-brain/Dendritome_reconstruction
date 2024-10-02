//
// Created by muyezhu on 10/13/17.
//

#ifndef MCP3D_IMAGE_UTILS_HPP
#define MCP3D_IMAGE_UTILS_HPP

#include <cstdint>
#include <type_traits>
#include <utility>
#include <array>
#include <stack>
#include <algorithm>
#include <opencv2/core/core.hpp>
#include "common/mcp3d_common.hpp"
#include "image_interface/mcp3d_voxel_types.hpp"
#include "image_interface/mcp3d_tiff_utils.hpp"
#include "mcp3d_image_view.hpp"


namespace mcp3d
{
/// d0 is the slowest varying dimension, d2 is the fastest. the *Robust version
/// performs more range checks. d0 - d2 are coordinates along each dimension.
/// dims_ contain size of each dimensions
int64_t LinearAddressRobust(const std::vector<int> &dims_, int d0, int d1, int d2 = 0);

int64_t LinearAddress(const std::vector<int> &dims_, int d0, int d1, int d2 = 0);

/// given the dimensions of a 3d volume and a linear address within the volume,
/// return zyx coordinates. address = z * ydim * xdim + y * xdim + xdim
std::array<int, 3> ZyxFromLinearAddress(const std::vector<int> &dims, int64_t address);

void ZyxFromLinearAddress(const std::vector<int> &dims, int64_t address, int* zyx);

/// connectivity_type = 1, 2, 3 respectively for 6, 18 and 26 neighbor connectivity
/// neighbor_addresses will be set to length 27 array
/// first entry of the array is the number of neighbor addresses. array positions
/// exceeding number of neighbor addresses will occupy are filled with -1
void NeighborAddresses(const std::vector<int> &dims, int64_t address, int connectivity_type,
                       std::unique_ptr<int64_t []>& neighbor_addresses);

void NeighborAddresses(const std::vector<int> &dims, int z, int y, int x,
                       int connectivity_type, std::unique_ptr<int64_t []>& neighbor_addresses);

/// number of strided steps within the given extent
/// e.g., extent = 5, stride = 1 => strided extent = 5
///       extent = 5, stride = 2 => strided extent = 3
template <typename T1 = int, typename T2 = int>
std::vector<int> StridedExtents(const std::vector<T1> &extents_, const std::vector<T2> &strides_ = std::vector<T2>{});

inline int StridedExtent(int extent, int stride)
{ return mcp3d::StridedExtents(std::vector<int>({extent}), std::vector<int>({stride}))[0]; }

/// number of voxels within a stride which is severed by boundaries of extents
/// for example, along a given dimension, extent = 5, stride = 2, gives a
/// residual stride of 1 (0 1 | 2 3 | 4 [boundary] 5)
template <typename T1 = int, typename T2 = int>
std::vector<int> ResidualStrides(const std::vector<T1> &extents_, const std::vector<T2> &strides_ = std::vector<T2>{});

int ResidualStride(int extent, int stride);

template <typename VType>
void SetConstant(void *data, const std::vector<int> &dims, VType value = VType{0});

template <typename VType>
void SetRandom(void *data, const std::vector<int> &dims);

template <typename VType>
void SetConstantBlock(void *data, const std::vector<int> &dims, const MImageBlock &image_block_ = mcp3d::MImageBlock{}, VType value = VType{0});

template <typename VType>
void SetConstantMultiBlocks(void *data, const std::vector<int> &dims, const std::vector<mcp3d::MImageBlock> &image_blocks, VType value = VType{0});

/// copy volumes of data from src to dst. the desired behavior is to copy
/// src[x_offset: x_offset + x_extent: x_stride,
///     y_offset: y_offset + y_extent: y_stride,
///     z_offset: z_offset + z_extent: z_stride] to
/// dst as continugous memory
/// calls CopyDataVolume(void*, const std::vector<int> &, void*, const std::vector<int> &,
///                      const mcp3d::MImageBlock&, const mcp3d::MImageBlock&)
/// will allocate necessary memory for dst if dst is null
template<typename VType>
void CopyDataVolume(void* dst, void *src, const std::vector<int> &src_dims, const MImageBlock& src_block = MImageBlock{});

/// copy volumes of data from src to dst. the desired behavior is to copy
/// src[x_offset: x_offset + x_extent: x_stride,
///     y_offset: y_offset + y_extent: y_stride,
///     z_offset: z_offset + z_extent: z_stride] to
/// dst[x_offset: x_offset + x_extent: x_stride,
///     y_offset: y_offset + y_extent: y_stride,
///     z_offset: z_offset + z_extent: z_stride]
template<typename VType>
void CopyDataVolume(void *dst, const std::vector<int> &dst_dims, void *src, const std::vector<int> &src_dims,
                    const MImageBlock& dst_block = MImageBlock{}, const MImageBlock& src_block = MImageBlock{});

/// does not require templated dst and src type, but need size (number of bytes)
/// per element in the array layout. useful when type size is known
void CopyDataVolume(void *dst, const std::vector<int> &dst_dims, void *src, const std::vector<int> &src_dims, int type_size,
                    const MImageBlock& dst_block = MImageBlock{}, const MImageBlock& src_block = MImageBlock{});

/// the two pointers must have the same type and same underlying array dimensions
bool VolumeEqual(const void *ptr1, const void *ptr2, const std::vector<int> &dims, int element_size);

/// array implementation of quicksort. sorting is in place
template <typename VType>
void QuickSort(void* data, long n_elements);

}

template <typename T1, typename T2>
std::vector<int> mcp3d::StridedExtents(const std::vector<T1> &extents_, const std::vector<T2> &strides_)
{
    static_assert(std::is_integral<T1>() && std::is_integral<T2>(),
                  "extents and strides must be integral types");
    MCP3D_ASSERT(!extents_.empty())
    MCP3D_ASSERT(mcp3d::AllPositive(extents_))
    std::vector<T1> extents = extents_;
    std::vector<T2> strides = strides_;
    if (strides.empty())
        strides.insert(strides.end(), extents.size(), T2{1});
    MCP3D_ASSERT(mcp3d::AllPositive(strides))
    MCP3D_ASSERT(extents.size() == strides.size())
    std::vector<int> strided_extents;
    for (size_t i = 0; i < extents.size(); ++i)
        strided_extents.push_back(extents[i] / strides[i] + (int)(extents[i] % strides[i] > 0));
    return strided_extents;
};

template <typename T1, typename T2>
std::vector<int> mcp3d::ResidualStrides(const std::vector<T1> &extents_, const std::vector<T2> &strides_)
{
    static_assert(std::is_integral<T1>() && std::is_integral<T2>(),
                  "extents and strides must be integral types");
    MCP3D_ASSERT(!extents_.empty())
    std::vector<T1> extents = extents_;
    std::vector<T2> strides = strides_;
    if (strides.empty())
        strides.insert(strides.end(), extents.size(), T2{1});
    MCP3D_ASSERT(extents.size() == strides.size())
    if (sizeof(T1) > sizeof(int) || sizeof(T2) > sizeof(int))
        MCP3D_MESSAGE("warning: input types may not fit in int type")
    std::vector<int> residual;
    for (size_t i = 0; i < extents.size(); ++i)
        residual.push_back(mcp3d::ResidualStride((int)extents[i], (int)strides[i]));
    return residual;
}

template <typename VType>
void mcp3d::SetConstant(void *data, const std::vector<int> &dims, VType value)
{
    MCP3D_ASSERT(data)
    static_assert(std::is_arithmetic<VType>(), "must have arithmetic element types");
    MCP3D_ASSERT(!dims.empty())
    MCP3D_ASSERT(mcp3d::AllPositive(dims))
    if (value == (VType)0)
        memset(data, value, mcp3d::ReduceProdSeq<size_t>(dims) * sizeof(VType));
    else
    {
        long addr_end = mcp3d::ReduceProdSeq<long>(dims);
        for (long i = 0; i < addr_end; ++i)
            ((VType*)data)[i] = value;
    }
}

template <typename VType>
void mcp3d::SetRandom(void *data, const std::vector<int> &dims)
{
    MCP3D_ASSERT(data)
    static_assert(std::is_arithmetic<VType>(), "must have arithmetic element types");
    MCP3D_ASSERT(mcp3d::AllPositive(dims))
    auto n = mcp3d::ReduceProdSeq<size_t>(dims);
    MCPTensor1DMap<VType> map((VType*)data, n);
    map.setRandom();
}

template <typename VType>
void mcp3d::SetConstantBlock(void *data, const std::vector<int> &dims, const mcp3d::MImageBlock &image_block_, VType value)
{
    MCP3D_ASSERT(data)
    static_assert(std::is_arithmetic<VType>(),
                  "must have arithmetic element types");
    MCP3D_ASSERT(dims.size() == 3)
    if (image_block_.empty())
    {
        mcp3d::SetConstant<VType>(data, dims, value);
        return;
    }
    std::vector<int> offsets(move(image_block_.offsets())),
                     extents(move(image_block_.extents())),
                     strides(move(image_block_.strides()));
    for (int i = 0; i < 3; ++i)
        if (extents[i] == 0)
            extents[i] = dims[i];

    MCP3D_ASSERT(AllNonNegative(offsets) && AllNonNegative(extents))
    MCP3D_ASSERT(AllPositive(strides))
    int z_end = std::min(dims[0], offsets[0] + extents[0]);
    int y_end = std::min(dims[1], offsets[1] + extents[1]);
    int x_end = std::min(dims[2], offsets[2] + extents[2]);
    if (offsets == std::vector<int>({0, 0, 0}) &&
        std::vector<int>({y_end, y_end, x_end}) == dims)
        mcp3d::SetConstant(data, dims, value);
    for (int z = offsets[0]; z < z_end; z += strides[0])
    {
        long addr;
        if (value == 0 && mcp3d::ReduceProdSeq<int>(strides) == 1)
        {
            size_t n_bytes = sizeof(VType) * extents[2];
            for (int y = offsets[1]; y < y_end; y += strides[1])
            {
                addr = mcp3d::LinearAddress(dims, z, offsets[1], offsets[0]);
                memset(((VType*)data) + addr, 0, n_bytes);
            }
        }
        else
        {
            for (int y = offsets[1]; y < y_end; y += strides[1])
            {
                int x = offsets[2];
                addr = mcp3d::LinearAddress(dims, z, y, x);
                for (; x < x_end; x += strides[2])
                    ((VType*)data)[addr++] = value;
            }
        }
    }
}

template <typename VType>
void mcp3d::SetConstantMultiBlocks(void *data, const std::vector<int> &dims, const std::vector<mcp3d::MImageBlock> &image_blocks, VType value)
{
    MCP3D_ASSERT(data)
    static_assert(std::is_arithmetic<VType>(),
                  "must have arithmetic element types");
    if (image_blocks.empty())
        SetConstant<VType>(data, dims, value);
    for (const mcp3d::MImageBlock& block: image_blocks)
        if (block.empty())
        {
            SetConstant<VType>(data, dims, value);
            return;
        }
    for (const mcp3d::MImageBlock& block: image_blocks)
        SetConstantBlock<VType>(data, dims, block, value);
}


template <typename VType>
void mcp3d::CopyDataVolume(void* dst, void *src, const std::vector<int> &src_dims, const MImageBlock &src_block)
{
    static_assert(std::is_arithmetic<VType>(), "must have arithmetic element types");
    MCP3D_ASSERT(src)
    std::vector<int> src_extents(move(src_block.extents())),
                     src_strides(move(src_block.strides()));
    // for any missing extents value (0 value), fill it in as if
    // offsets: end of dimension was given
    for (int i = 0; i < 3; ++i)
        if (src_extents[i] == 0)
            src_extents[i] = src_dims[i] - src_block.offsets()[i];
    std::vector<int> dst_dims = mcp3d::StridedExtents(src_extents, src_strides);
    auto vtype_dst = (VType*) dst;
    if (!vtype_dst)
        vtype_dst = new (std::nothrow) VType[mcp3d::ReduceProdSeq<size_t>(dst_dims)];
    mcp3d::MImageBlock dst_block({0, 0, 0}, dst_dims, {1, 1, 1});
    mcp3d::CopyDataVolume<VType>(vtype_dst, dst_dims, src, src_dims, dst_block, src_block);
}

template<typename VType>
void mcp3d::CopyDataVolume(void *dst, const std::vector<int> &dst_dims, void *src, const std::vector<int> &src_dims,
                           const mcp3d::MImageBlock& dst_block, const mcp3d::MImageBlock& src_block)
{
    static_assert(std::is_arithmetic<VType>(), "must have arithmetic element types");
    MCP3D_ASSERT(dst && src)
    MCP3D_ASSERT(dst_dims.size() == 3 && src_dims.size() == 3)
    MCP3D_ASSERT(AllPositive(dst_dims) && AllPositive(src_dims))
    std::vector<int> src_offsets(move(src_block.offsets())),
                     src_extents(move(src_block.extents())),
                     src_strides(move(src_block.strides())),
                     dst_offsets(move(dst_block.offsets())),
                     dst_strides(move(dst_block.strides()));
    // for any missing extents value (0 value), fill it in as if
    // offsets: end of dimension was given
    for (int i = 0; i < 3; ++i)
        if (src_extents[i] == 0)
            src_extents[i] = src_dims[i] - src_offsets[i];

    MCP3D_ASSERT(AllNonNegative(src_offsets) && AllNonNegative(src_extents))
    MCP3D_ASSERT(AllPositive(src_strides))
    MCP3D_ASSERT(AllNonNegative(SubtractSeq<int>(src_dims, AddSeq<int>(src_offsets, src_extents))))
    if (ReduceProdSeq<int>(src_extents) == 0)
        return;

    MCP3D_ASSERT(AllNonNegative(dst_offsets))
    MCP3D_ASSERT(AllPositive(dst_strides))
    std::vector<int> src_extents_strided =
            mcp3d::StridedExtents(src_extents, src_strides);
    std::vector<int> dst_extents = mcp3d::MultiplySeq<int>(dst_strides,
                                                           src_extents_strided);
    std::vector<int> dst_remain_extents = SubtractSeq<int>(dst_dims, AddSeq<int>(dst_offsets, dst_extents));
    MCP3D_ASSERT(AllNonNegative(dst_remain_extents))

    // copy data
    for (int z_dst = dst_offsets[0], z_src = src_offsets[0];
         z_src < src_offsets[0] + src_extents[0];
         z_dst += dst_strides[0], z_src += src_strides[0])
    {
        // memcpy for simple case
        if (mcp3d::ReduceProdSeq<int>(src_strides) == 1 && mcp3d::ReduceProdSeq<int>(dst_strides) == 1)
        {
            if (mcp3d::AllZeros(src_offsets) && mcp3d::AllZeros(dst_offsets) &&
                src_extents == src_dims && dst_extents == dst_dims && src_dims == dst_dims)
            {
                #ifdef VERBOSE
                if (z_dst == dst_offsets[2])
                    MCP3D_MESSAGE("copy continuous block to same type")
                #endif
                size_t n_bytes = sizeof(VType) * mcp3d::ReduceProdSeq<int>(src_dims);
                memcpy(dst, src, n_bytes);
            }
            else
            {
                #ifdef VERBOSE
                if (z_dst == dst_offsets[2])
                    MCP3D_MESSAGE("copy discontinuous unit stride block to sampe type");
                #endif
                long addr_dst = mcp3d::LinearAddress(dst_dims, z_dst, dst_offsets[1], dst_offsets[2]),
                     addr_src = mcp3d::LinearAddress(src_dims, z_src, src_offsets[1], src_offsets[2]);
                for(int y_src = src_offsets[1]; y_src < src_offsets[1] + src_extents[1]; ++y_src)
                {
                    size_t n_bytes = src_extents[2] * sizeof(VType);
                    memcpy(((VType *) dst) + addr_dst,
                           ((VType *) src) + addr_src, n_bytes);
                    addr_dst += dst_dims[2];
                    addr_src += src_dims[2];
                }
            }
        }
        else
        {
            #ifdef VERBOSE
            if (z_dst == dst_offsets[2])
                MCP3D_MESSAGE("general copy")
            #endif
            for (int y_dst = dst_offsets[1], y_src = src_offsets[1];
                 y_src < src_offsets[1] + src_extents[1];
                 y_dst += dst_strides[1], y_src += src_strides[1])
            {
                int x_dst = dst_offsets[2], x_src = src_offsets[2];
                long addr_src = mcp3d::LinearAddress(src_dims, z_src, y_src, x_src),
                     addr_dst = mcp3d::LinearAddress(dst_dims, z_dst, y_dst, x_dst);
                for (; x_src < src_offsets[2] + src_extents[2]; x_src += src_strides[2])
                {
                    (reinterpret_cast<VType*>(dst))[addr_dst] = (reinterpret_cast<VType*>(src))[addr_src];
                    addr_dst += dst_strides[2];
                    addr_src += src_strides[2];
                }
            }
        }
    }
}

template <typename VType>
long QuickSortPartition(void* data, long low, long high)
{
    static_assert(std::is_arithmetic<VType>(), "must have arithmetic element types");
    auto vtype_data = (VType*)data;
    VType pivot = vtype_data[low + (high - low) / 2];
    while (true)
    {
        while (vtype_data[low] < pivot)
            ++low;
        while (vtype_data[high] > pivot)
            --high;
        if (low >= high)
            return high;
        std::swap<VType>(vtype_data[low], vtype_data[high]);
        ++low;
        --high;
    }
}

template <typename VType>
void QuickSort(void *data, long low, long high)
{
    std::stack<long> index_stack;
    index_stack.push(high);
    index_stack.push(low);
    long low_index, high_index, pivot_index;
    while (!index_stack.empty())
    {
        low_index = index_stack.top();
        index_stack.pop();
        high_index = index_stack.top();
        index_stack.pop();
        if (low_index < high_index)
        {
            pivot_index = ::QuickSortPartition<VType>(data, low_index, high_index);
            index_stack.push(pivot_index);
            index_stack.push(low_index);
            index_stack.push(high_index);
            index_stack.push(pivot_index + 1);
        }
    }
}

template <typename VType>
void mcp3d::QuickSort(void *data, long n_elements)
{
    if (n_elements <= 1)
        return;
    ::QuickSort<VType>(data, 0, n_elements - 1);
}


#endif //MCP3D_IMAGE_UTILS_HPP
